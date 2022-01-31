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
    data-toggle="buttons-radio"---
title: "Fetching JSON data from REST APIs"
date: "2015-09-06"
output:
  html_document
vignette: >
  %\VignetteIndexEntry{Fetching JSON data from REST APIs}
  %\VignetteEngine{knitr::rmarkdown}
  \usepackage[utf8]{inputenc}
---



This section lists some examples of public HTTP APIs that publish data in JSON format. These are great to get a sense of the complex structures that are encountered in real world JSON data. All services are free, but some require registration/authentication. Each example returns lots of data, therefore not all output is printed in this document.


```r
library(jsonlite)
```

## Github

Github is an online code repository and has APIs to get live data on almost all activity. Below some examples from a well known R package and author:


```r
hadley_orgs <- fromJSON("https://api.github.com/users/hadley/orgs")
hadley_repos <- fromJSON("https://api.github.com/users/hadley/repos")
gg_commits <- fromJSON("https://api.github.com/repos/hadley/ggplot2/commits")
gg_issues <- fromJSON("https://api.github.com/repos/hadley/ggplot2/issues")

#latest issues
paste(format(gg_issues$user$login), ":", gg_issues$title)
```

```
 [1] "idavydov     : annotate(\"segment\") wrong position if limits are inverted"                      
 [2] "ben519       : geom_polygon doesn't make NA values grey when using continuous fill"              
 [3] "has2k1       : Fix multiple tiny issues in the position classes"                                 
 [4] "neggert      : Problem with geom_bar position=fill and faceting"                                 
 [5] "robertzk     : Fix typo in geom_linerange docs."                                                 
 [6] "lionel-      : stat_bar() gets confused with numeric discrete data?"                             
 [7] "daattali     : Request: support theme axis.ticks.length.x and axis.ticks.length.y"               
 [8] "sethchandler : Documentation error on %+replace% ?"                                              
 [9] "daattali     : dev version 1.0.1.9003 has some breaking changes"                                 
[10] "lionel-      : Labels"                                                                           
[11] "nutterb      : legend for `geom_line` colour disappears when `alpha` < 1.0"                      
[12] "wch          : scale_name property should be removed from Scale objects"                         
[13] "wch          : scale_details arguments in Coords should be renamed panel_scales or scale"        
[14] "wch          : ScalesList-related functions should be moved into ggproto object"                 
[15] "wch          : update_geom_defaults and update_stat_defaults should accept Geom and Stat objects"
[16] "wch          : Make some ggproto objects immutable. Closes #1237"                                
[17] "and3k        : Control size of the border and padding of geom_label"                             
[18] "hadley       : Consistent argument order and formatting for layer functions"                     
[19] "hadley       : Consistently handle missing values"                                               
[20] "cmohamma     : fortify causes fatal error"                                                       
[21] "lionel-      : Flawed `label_bquote()` implementation"                                           
[22] "beroe        : Create alias for `colors=` in `scale_color_gradientn()`"                          
[23] "and3k        : hjust broken in y facets"                                                         
[24] "joranE       : Allow color bar guides for alpha scales"                                          
[25] "hadley       : dir = \"v\" also needs to swap nrow and ncol"                                     
[26] "joranE       : Add examples for removing guides"                                                 
[27] "lionel-      : New approach for horizontal layers"                                               
[28] "bbolker      : add horizontal linerange geom"                                                    
[29] "hadley       : Write vignette about grid"                                                        
[30] "hadley       : Immutable flag for ggproto objects"                                               
```

## CitiBike NYC

A single public API that shows location, status and current availability for all stations in the New York City bike sharing imitative.


```r
citibike <- fromJSON("http://citibikenyc.com/stations/json")
stations <- citibike$stationBeanList
colnames(stations)
```

```
 [1] "id"                    "stationName"          
 [3] "availableDocks"        "totalDocks"           
 [5] "latitude"              "longitude"            
 [7] "statusValue"           "statusKey"            
 [9] "availableBikes"        "stAddress1"           
[11] "stAddress2"            "city"                 
[13] "postalCode"            "location"             
[15] "altitude"              "testStation"          
[17] "lastCommunicationTime" "landMark"             
```

```r
nrow(stations)
```

```
[1] 509
```

## Ergast

The Ergast Developer API is an experimental web service which provides a historical record of motor racing data for non-commercial purposes.


```r
res <- fromJSON('http://ergast.com/api/f1/2004/1/results.json')
drivers <- res$MRData$RaceTable$Races$Results[[1]]$Driver
colnames(drivers)
```

```
[1] "driverId"        "code"            "url"             "givenName"      
[5] "familyName"      "dateOfBirth"     "nationality"     "permanentNumber"
```

```r
drivers[1:10, c("givenName", "familyName", "code", "nationality")]
```

```
   givenName    familyName code nationality
1    Michael    Schumacher  MSC      German
2     Rubens   Barrichello  BAR   Brazilian
3   Fernando        Alonso  ALO     Spanish
4       Ralf    Schumacher  SCH      German
5       Juan Pablo Montoya  MON   Colombian
6     Jenson        Button  BUT     British
7      Jarno        Trulli  TRU     Italian
8      David     Coulthard  COU     British
9     Takuma          Sato  SAT    Japanese
10 Giancarlo    Fisichella  FIS     Italian
```


## ProPublica

Below an example from the [ProPublica Nonprofit Explorer API](http://projects.propublica.org/nonprofits/api) where we retrieve the first 10 pages of tax-exempt organizations in the USA, ordered by revenue. The `rbind.pages` function is used to combine the pages into a single data frame.



```r
#store all pages in a list first
baseurl <- "https://projects.propublica.org/nonprofits/api/v1/search.json?order=revenue&sort_order=desc"
pages <- list()
for(i in 0:10){
  mydata <- fromJSON(paste0(baseurl, "&page=", i), flatten=TRUE)
  message("Retrieving page ", i)
  pages[[i+1]] <- mydata$filings
}

#combine all into one
filings <- rbind.pages(pages)

#check output
nrow(filings)
```

```
[1] 275
```

```r
filings[1:10, c("organization.sub_name", "organization.city", "totrevenue")]
```

```
                              organization.sub_name organization.city
1                 KAISER FOUNDATION HEALTH PLAN INC           OAKLAND
2                 KAISER FOUNDATION HEALTH PLAN INC           OAKLAND
3                 KAISER FOUNDATION HEALTH PLAN INC           OAKLAND
4  DAVIDSON COUNTY COMMUNITY COLLEGE FOUNDATION INC         LEXINGTON
5                       KAISER FOUNDATION HOSPITALS           OAKLAND
6                       KAISER FOUNDATION HOSPITALS           OAKLAND
7                       KAISER FOUNDATION HOSPITALS           OAKLAND
8                   PARTNERS HEALTHCARE SYSTEM INC        CHARLESTOWN
9                   PARTNERS HEALTHCARE SYSTEM INC        CHARLESTOWN
10                  PARTNERS HEALTHCARE SYSTEM INC        CHARLESTOWN
    totrevenue
1  42346486950
2  40148558254
3  37786011714
4  30821445312
5  20013171194
6  18543043972
7  17980030355
8  10619215354
9  10452560305
10  9636630380
```


## New York Times

The New York Times has several APIs as part of the NYT developer network. These interface to data from various departments, such as news articles, book reviews, real estate, etc. Registration is required (but free) and a key can be obtained at [here](http://developer.nytimes.com/docs/reference/keys). The code below includes some example keys for illustration purposes.


```r
#search for articles
article_key <- "&api-key=c2fede7bd9aea57c898f538e5ec0a1ee:6:68700045"
url <- "http://api.nytimes.com/svc/search/v2/articlesearch.json?q=obamacare+socialism"
req <- fromJSON(paste0(url, article_key))
articles <- req$response$docs
colnames(articles)
```

```
 [1] "web_url"          "snippet"          "lead_paragraph"  
 [4] "abstract"         "print_page"       "blog"            
 [7] "source"           "multimedia"       "headline"        
[10] "keywords"         "pub_date"         "document_type"   
[13] "news_desk"        "section_name"     "subsection_name" 
[16] "byline"           "type_of_material" "_id"             
[19] "word_count"      
```

```r
#search for best sellers
bestseller_key <- "&api-key=5e260a86a6301f55546c83a47d139b0d:3:68700045"
url <- "http://api.nytimes.com/svc/books/v2/lists/overview.json?published_date=2013-01-01"
req <- fromJSON(paste0(url, bestseller_key))
bestsellers <- req$results$list
category1 <- bestsellers[[1, "books"]]
subset(category1, select = c("author", "title", "publisher"))
```

```
           author                title                  publisher
1   Gillian Flynn            GONE GIRL           Crown Publishing
2    John Grisham        THE RACKETEER Knopf Doubleday Publishing
3       E L James FIFTY SHADES OF GREY Knopf Doubleday Publishing
4 Nicholas Sparks           SAFE HAVEN   Grand Central Publishing
5  David Baldacci        THE FORGOTTEN   Grand Central Publishing
```

```r
#movie reviews
movie_key <- "&api-key=5a3daaeee6bbc6b9df16284bc575e5ba:0:68700045"
url <- "http://api.nytimes.com/svc/movies/v2/reviews/dvd-picks.json?order=by-date"
req <- fromJSON(paste0(url, movie_key))
reviews <- req$results
colnames(reviews)
```

```
 [1] "nyt_movie_id"     "display_title"    "sort_name"       
 [4] "mpaa_rating"      "critics_pick"     "thousand_best"   
 [7] "byline"           "headline"         "capsule_review"  
[10] "summary_short"    "publication_date" "opening_date"    
[13] "dvd_release_date" "date_updated"     "seo_name"        
[16] "link"             "related_urls"     "multimedia"      
```

```r
reviews[1:5, c("display_title", "byline", "mpaa_rating")]
```

```
       display_title         byline mpaa_rating
1    Tom at the Farm Stephen Holden          NR
2     A Little Chaos Stephen Holden           R
3           Big Game   Andy Webster        PG13
4          Balls Out   Andy Webster           R
5 Mad Max: Fury Road    A. O. Scott           R
```

## CrunchBase

CrunchBase is the free database of technology companies, people, and investors that anyone can edit.


```r
key <- "f6dv6cas5vw7arn5b9d7mdm3"
res <- fromJSON(paste0("http://api.crunchbase.com/v/1/search.js?query=R&api_key=", key))
head(res$results)
```

## Sunlight Foundation

The Sunlight Foundation is a non-profit that helps to make government transparent and accountable through data, tools, policy and journalism. Register a free key at [here](http://sunlightfoundation.com/api/accounts/register/). An example key is provided.


```r
key <- "&apikey=39c83d5a4acc42be993ee637e2e4ba3d"

#Find bills about drones
drone_bills <- fromJSON(paste0("http://openstates.org/api/v1/bills/?q=drone", key))
drone_bills$title <- substring(drone_bills$title, 1, 40)
print(drone_bills[1:5, c("title", "state", "chamber", "type")])
```

```
                                     title state chamber type
1                            WILDLIFE-TECH    il   lower bill
2 Criminalizes the unlawful use of an unma    ny   lower bill
3 Criminalizes the unlawful use of an unma    ny   lower bill
4 Relating to: criminal procedure and prov    wi   lower bill
5 Relating to: criminal procedure and prov    wi   upper bill
```

```r
#Congress mentioning "constitution"
res <- fromJSON(paste0("http://capitolwords.org/api/1/dates.json?phrase=immigration", key))
wordcount <- res$results
wordcount$day <- as.Date(wordcount$day)
summary(wordcount)
```

```
     count              day               raw_count      
 Min.   :   1.00   Min.   :1996-01-02   Min.   :   1.00  
 1st Qu.:   3.00   1st Qu.:2001-01-22   1st Qu.:   3.00  
 Median :   8.00   Median :2005-11-16   Median :   8.00  
 Mean   :  25.27   Mean   :2005-10-02   Mean   :  25.27  
 3rd Qu.:  21.00   3rd Qu.:2010-05-12   3rd Qu.:  21.00  
 Max.   :1835.00   Max.   :2015-08-05   Max.   :1835.00  
```

```r
#Local legislators
legislators <- fromJSON(paste0("http://congress.api.sunlightfoundation.com/",
  "legislators/locate?latitude=42.96&longitude=-108.09", key))
subset(legislators$results, select=c("last_name", "chamber", "term_start", "twitter_id"))
```

```
  last_name chamber term_start      twitter_id
1    Lummis   house 2015-01-06   CynthiaLummis
2      Enzi  senate 2015-01-06     SenatorEnzi
3  Barrasso  senate 2013-01-03 SenJohnBarrasso
```

## Twitter

The twitter API requires OAuth2 authentication. Some example code:


```r
#Create your own appication key at https://dev.twitter.com/apps
consumer_key = "EZRy5JzOH2QQmVAe9B4j2w";
consumer_secret = "OIDC4MdfZJ82nbwpZfoUO4WOLTYjoRhpHRAWj6JMec";

#Use basic auth
library(httr)
secret <- RCurl::base64(paste(consumer_key, consumer_secret, sep = ":"));
req <- POST("https://api.twitter.com/oauth2/token",
  add_headers(
    "Authorization" = paste("Basic", secret),
    "Content-Type" = "application/x-www-form-urlencoded;charset=UTF-8"
  ),
  body = "grant_type=client_credentials"
);

#Extract the access token
token <- paste("Bearer", content(req)$access_token)

#Actual API call
url <- "https://api.twitter.com/1.1/statuses/user_timeline.json?count=10&screen_name=Rbloggers"
req <- GET(url, add_headers(Authorization = token))
json <- content(req, as = "text")
tweets <- fromJSON(json)
substring(tweets$text, 1, 100)
```

```
 [1] "Analysing longitudinal data: Multilevel growth models (II) http://t.co/unUxszG7VJ #rstats"           
 [2] "RcppDE 0.1.4 http://t.co/3qPhFzoOpj #rstats"                                                         
 [3] "Minimalist Maps http://t.co/fpkNznuCoX #rstats"                                                      
 [4] "Tutorials freely available of course I taught: including ggplot2, dplyr and shiny http://t.co/WsxX4U"
 [5] "Deploying Shiny apps with shinyapps.io http://t.co/tjef1pbKLt #rstats"                               
 [6] "Bootstrap Evaluation of Clusters http://t.co/EbY7ziKCz5 #rstats"                                     
 [7] "Add external code to Rmarkdown http://t.co/RCJEmS8gyP #rstats"                                       
 [8] "Linear models with weighted observations http://t.co/pUoHpvxAGC #rstats"                             
 [9] "dplyr 0.4.3 http://t.co/ze3zc8t7qj #rstats"                                                          
[10] "xkcd survey and the power to shape the internet http://t.co/vNaKhxWxE4 #rstats"                      
```

---
Title: "Getting started with JSON and jsonlite"
date: "`r Sys.Date()`"
output:
  html_document
vignette: >
  %\VignetteIndexEntry{Getting started with JSON and jsonlite}
  %\VignetteEngine{knitr::rmarkdown}
  \usepackage[utf8]{inputenc}
---


```{r echo=FALSE}
library(knitr)
opts_chunk$set(comment="")

#this replaces tabs by spaces because latex-verbatim doesn't like tabs
#no longer needed because yajl does not use tabs.
#toJSON <- function(...){
#  gsub("\t", "  ", jsonlite::toJSON(...), fixed=TRUE);
#}
```

# Getting started with JSON and jsonlite

The jsonlite package is a JSON parser/generator optimized for the web. Its main strength is that it implements a bidirectional mapping between JSON data and the most important R data types. Thereby we can convert between R objects and JSON without loss of type or information, and without the need for any manual data munging. This is ideal for interacting with web APIs, or to build pipelines where data structures seamlessly flow in and out of R using JSON.

```{r message=FALSE}
library(jsonlite)
all.equal(mtcars, fromJSON(toJSON(mtcars)))
```

This vignette introduces basic concepts to get started with jsonlite. For a more detailed outline and motivation of the mapping, see: [arXiv:1403.2805](http://arxiv.org/abs/1403.2805).

## Simplification

Simplification is the process where JSON arrays automatically get converted from a list into a more specific R class. The `fromJSON` function has 3 arguments which control the simplification process: `simplifyVector`, `simplifyDataFrame` and `simplifyMatrix`. Each one is enabled by default.

| JSON structure        | Example JSON data                                        | Simplifies to R class | Argument in fromJSON | 
| ----------------------|----------------------------------------------------------|-----------------------|----------------------|
| Array of primitives   | `["Amsterdam", "Rotterdam", "Utrecht", "Den Haag"]`      | Atomic Vector         | simplifyVector       | 
| Array of objects      | `[{"name":"Erik", "age":43}, {"name":"Anna", "age":32}]` | Data Frame            | simplifyDataFrame    | 
| Array of arrays       | `[ [1, 2, 3], [4, 5, 6] ]`                               | Matrix                | simplifyMatrix       |

### Atomic Vectors

When `simplifyVector` is enabled, JSON arrays containing **primitives** (strings, numbers, booleans or null) simplify into an atomic vector:

```{r}
# A JSON array of primitives
json <- '["Mario", "Peach", null, "Bowser"]'

# Simplifies into an atomic vector
fromJSON(json)
```

Without simplification, any JSON array turns into a list: 

```{r}
# No simplification:
fromJSON(json, simplifyVector = FALSE)
```


### Data Frames

When `simplifyDataFrame` is enabled, JSON arrays containing **objects** (key-value pairs) simplify into a data frame:

```{r}
json <-
'[
  {"Name" : "Mario", "Age" : 32, "Occupation" : "Plumber"}, 
  {"Name" : "Peach", "Age" : 21, "Occupation" : "Princess"},
  {},
  {"Name" : "Bowser", "Occupation" : "Koopa"}
]'
mydf <- fromJSON(json)
mydf
```

The data frame gets converted back into the original JSON structure by `toJSON` (whitespace and line breaks are ignorable in JSON).

```{r}
mydf$Ranking <- c(3, 1, 2, 4)
toJSON(mydf, pretty=TRUE)
```

Hence you can go back and forth between dataframes and JSON, without any manual data restructuring.

### Matrices and Arrays

When `simplifyMatrix` is enabled, JSON arrays containing **equal-length sub-arrays** simplify into a matrix (or higher order R array):

```{r}
json <- '[
  [1, 2, 3, 4],
  [5, 6, 7, 8],
  [9, 10, 11, 12]
]'
mymatrix <- fromJSON(json)
mymatrix
```

Again, we can use `toJSON` to convert the matrix or array back into the original JSON structure:

```{r}
toJSON(mymatrix, pretty = TRUE)
```

The simplification works for arrays of arbitrary dimensionality, as long as the dimensions match (R does not support ragged arrays).

```{r}
json <- '[
   [[1, 2], 
    [3, 4]],
   [[5, 6], 
    [7, 8]],
   [[9, 10],
    [11, 12]]
]'
myarray <- fromJSON(json)
myarray[1, , ]
myarray[ , ,1]
```

This is all there is to it! For a more detailed outline and motivation of the mapping, see: [arXiv:1403.2805](http://arxiv.org/abs/1403.2805).
---
title: "Combining pages of JSON data with jsonlite"
date: "2015-09-06"
output:
  html_document
vignette: >
  %\VignetteIndexEntry{Combining pages of JSON data with jsonlite}
  %\VignetteEngine{knitr::rmarkdown}
  \usepackage[utf8]{inputenc}
---






The [jsonlite](https://cran.r-project.org/package=jsonlite) package is a `JSON` parser/generator for R which is optimized for pipelines and web APIs. It is used by the OpenCPU system and many other packages to get data in and out of R using the `JSON` format.

## A bidirectional mapping

One of the main strengths of `jsonlite` is that it implements a bidirectional [mapping](http://arxiv.org/abs/1403.2805) between JSON and data frames. Thereby it can convert nested collections of JSON records, as they often appear on the web, immediately into the appropriate R structure. For example to grab some data from ProPublica we can simply use:


```r
library(jsonlite)
mydata <- fromJSON("https://projects.propublica.org/forensics/geos.json", flatten = TRUE)
View(mydata)
```

The `mydata` object is a data frame which can be used directly for modeling or visualization, without the need for any further complicated data manipulation.

## Paging with jsonlite

A question that comes up frequently is how to combine pages of data. Most web APIs limit the amount of data that can be retrieved per request. If the client needs more data than what can fits in a single request, it needs to break down the data into multiple requests that each retrieve a fragment (page) of data, not unlike pages in a book. In practice this is often implemented using a `page` parameter in the API. Below an example from the [ProPublica Nonprofit Explorer API](http://projects.propublica.org/nonprofits/api) where we retrieve the first 3 pages of tax-exempt organizations in the USA, ordered by revenue:


```r
baseurl <- "https://projects.propublica.org/nonprofits/api/v1/search.json?order=revenue&sort_order=desc"
mydata0 <- fromJSON(paste0(baseurl, "&page=0"), flatten = TRUE)
mydata1 <- fromJSON(paste0(baseurl, "&page=1"), flatten = TRUE)
mydata2 <- fromJSON(paste0(baseurl, "&page=2"), flatten = TRUE)

#The actual data is in the filings element
mydata0$filings[1:10, c("organization.sub_name", "organization.city", "totrevenue")]
```

```
                              organization.sub_name organization.city
1                 KAISER FOUNDATION HEALTH PLAN INC           OAKLAND
2                 KAISER FOUNDATION HEALTH PLAN INC           OAKLAND
3                 KAISER FOUNDATION HEALTH PLAN INC           OAKLAND
4  DAVIDSON COUNTY COMMUNITY COLLEGE FOUNDATION INC         LEXINGTON
5                       KAISER FOUNDATION HOSPITALS           OAKLAND
6                       KAISER FOUNDATION HOSPITALS           OAKLAND
7                       KAISER FOUNDATION HOSPITALS           OAKLAND
8                   PARTNERS HEALTHCARE SYSTEM INC        CHARLESTOWN
9                   PARTNERS HEALTHCARE SYSTEM INC        CHARLESTOWN
10                  PARTNERS HEALTHCARE SYSTEM INC        CHARLESTOWN
    totrevenue
1  42346486950
2  40148558254
3  37786011714
4  30821445312
5  20013171194
6  18543043972
7  17980030355
8  10619215354
9  10452560305
10  9636630380
```

To analyze or visualize these data, we need to combine the pages into a single dataset. We can do this with the `rbind.pages` function. Note that in this example, the actual data is contained by the `filings` field:


```r
#Rows per data frame
nrow(mydata0$filings)
```

```
[1] 25
```

```r
#Combine data frames
filings <- rbind.pages(
  list(mydata0$filings, mydata1$filings, mydata2$filings)
)

#Total number of rows
nrow(filings)
```

```
[1] 75
```

## Automatically combining many pages

We can write a simple loop that automatically downloads and combines many pages. For example to retrieve the first 20 pages with non-profits from the example above:


```r
#store all pages in a list first
baseurl <- "https://projects.propublica.org/nonprofits/api/v1/search.json?order=revenue&sort_order=desc"
pages <- list()
for(i in 0:20){
  mydata <- fromJSON(paste0(baseurl, "&page=", i))
  message("Retrieving page ", i)
  pages[[i+1]] <- mydata$filings
}

#combine all into one
filings <- rbind.pages(pages)

#check output
nrow(filings)
```

```
[1] 525
```

```r
colnames(filings)
```

```
  [1] "tax_prd"               "tax_prd_yr"           
  [3] "formtype"              "pdf_url"              
  [5] "updated"               "totrevenue"           
  [7] "totfuncexpns"          "totassetsend"         
  [9] "totliabend"            "pct_compnsatncurrofcr"
 [11] "tax_pd"                "subseccd"             
 [13] "unrelbusinccd"         "initiationfees"       
 [15] "grsrcptspublicuse"     "grsincmembers"        
 [17] "grsincother"           "totcntrbgfts"         
 [19] "totprgmrevnue"         "invstmntinc"          
 [21] "txexmptbndsproceeds"   "royaltsinc"           
 [23] "grsrntsreal"           "grsrntsprsnl"         
 [25] "rntlexpnsreal"         "rntlexpnsprsnl"       
 [27] "rntlincreal"           "rntlincprsnl"         
 [29] "netrntlinc"            "grsalesecur"          
 [31] "grsalesothr"           "cstbasisecur"         
 [33] "cstbasisothr"          "gnlsecur"             
 [35] "gnlsothr"              "netgnls"              
 [37] "grsincfndrsng"         "lessdirfndrsng"       
 [39] "netincfndrsng"         "grsincgaming"         
 [41] "lessdirgaming"         "netincgaming"         
 [43] "grsalesinvent"         "lesscstofgoods"       
 [45] "netincsales"           "miscrevtot11e"        
 [47] "compnsatncurrofcr"     "othrsalwages"         
 [49] "payrolltx"             "profndraising"        
 [51] "txexmptbndsend"        "secrdmrtgsend"        
 [53] "unsecurednotesend"     "retainedearnend"      
 [55] "totnetassetend"        "nonpfrea"             
 [57] "gftgrntsrcvd170"       "txrevnuelevied170"    
 [59] "srvcsval170"           "grsinc170"            
 [61] "grsrcptsrelated170"    "totgftgrntrcvd509"    
 [63] "grsrcptsadmissn509"    "txrevnuelevied509"    
 [65] "srvcsval509"           "subtotsuppinc509"     
 [67] "totsupp509"            "ein"                  
 [69] "organization"          "eostatus"             
 [71] "tax_yr"                "operatingcd"          
 [73] "assetcdgen"            "transinccd"           
 [75] "subcd"                 "grscontrgifts"        
 [77] "intrstrvnue"           "dividndsamt"          
 [79] "totexcapgn"            "totexcapls"           
 [81] "grsprofitbus"          "otherincamt"          
 [83] "compofficers"          "contrpdpbks"          
 [85] "totrcptperbks"         "totexpnspbks"         
 [87] "excessrcpts"           "totexpnsexempt"       
 [89] "netinvstinc"           "totaxpyr"             
 [91] "adjnetinc"             "invstgovtoblig"       
 [93] "invstcorpstk"          "invstcorpbnd"         
 [95] "totinvstsec"           "fairmrktvalamt"       
 [97] "undistribincyr"        "cmpmininvstret"       
 [99] "sec4940notxcd"         "sec4940redtxcd"       
[101] "infleg"                "contractncd"          
[103] "claimstatcd"           "propexchcd"           
[105] "brwlndmnycd"           "furngoodscd"          
[107] "paidcmpncd"            "trnsothasstscd"       
[109] "agremkpaycd"           "undistrinccd"         
[111] "dirindirintcd"         "invstjexmptcd"        
[113] "propgndacd"            "excesshldcd"          
[115] "grntindivcd"           "nchrtygrntcd"         
[117] "nreligiouscd"          "grsrents"             
[119] "costsold"              "totrcptnetinc"        
[121] "trcptadjnetinc"        "topradmnexpnsa"       
[123] "topradmnexpnsb"        "topradmnexpnsd"       
[125] "totexpnsnetinc"        "totexpnsadjnet"       
[127] "othrcashamt"           "mrtgloans"            
[129] "othrinvstend"          "fairmrktvaleoy"       
[131] "mrtgnotespay"          "tfundnworth"          
[133] "invstexcisetx"         "sect511tx"            
[135] "subtitleatx"           "esttaxcr"             
[137] "txwithldsrc"           "txpaidf2758"          
[139] "erronbkupwthld"        "estpnlty"             
[141] "balduopt"              "crelamt"              
[143] "tfairmrktunuse"        "distribamt"           
[145] "adjnetinccola"         "adjnetinccolb"        
[147] "adjnetinccolc"         "adjnetinccold"        
[149] "adjnetinctot"          "qlfydistriba"         
[151] "qlfydistribb"          "qlfydistribc"         
[153] "qlfydistribd"          "qlfydistribtot"       
[155] "valassetscola"         "valassetscolb"        
[157] "valassetscolc"         "valassetscold"        
[159] "valassetstot"          "qlfyasseta"           
[161] "qlfyassetb"            "qlfyassetc"           
[163] "qlfyassetd"            "qlfyassettot"         
[165] "endwmntscola"          "endwmntscolb"         
[167] "endwmntscolc"          "endwmntscold"         
[169] "endwmntstot"           "totsuprtcola"         
[171] "totsuprtcolb"          "totsuprtcolc"         
[173] "totsuprtcold"          "totsuprttot"          
[175] "pubsuprtcola"          "pubsuprtcolb"         
[177] "pubsuprtcolc"          "pubsuprtcold"         
[179] "pubsuprttot"           "grsinvstinca"         
[181] "grsinvstincb"          "grsinvstincc"         
[183] "grsinvstincd"          "grsinvstinctot"       
```

From here, we can go straight to analyzing the filings data without any further tedious data manipulation.
\name{NEWS}
\title{News for Package 'rhdf5'}

\section{Changes in version 2.10.0}{
  \subsection{NEW FEATURES}{
    \itemize{
      \item Added support for HDF5 property lists.
      \item Added property list arguments to H5Dcreate and H5Dopen.
      \item New function h5readAttributes implemented that reads all HDF5 attributes of one object.
      \item New function h5version implemented.
      \item fillValue parameter added to h5createDataset.
      \item New low level general library functions H5Lcreate_external, H5Fis_hdf5, H5Fget_filesize, H5Fget_name, H5Pcreate, H5Pcopy, H5Pget_class, H5Pclose, H5Pclose_class, H5Pset_char_encoding, H5Pset_create_intermediate_group, H5Pset_chunk_cache, H5Pset_layout, H5Pset_chunk, H5Pget_chunk, H5Pset_deflate, H5Pset_fill_value, H5Pset_fill_time, H5Pset_alloc_time, H5Pequal implemented.
      \item Support for parallel Make (make -j)
    }
  }
  \subsection{USER VISIBLE CHANGES}{
    \itemize{
      \item A warning is shown in high level function (h5read, h5write and others), if an open HDF5 handle already exists for the specified filename.
    }
  }
  \subsection{BUG FIXES}{
    \itemize{
      \item Error in h5write for 0-length objects, as a consequence of automatic determining chunk size
      \item missing size parameter message in h5createDataset now correctly display
      \item checking for open file identifiers in h5read and h5ls now only searches for file names in open files, groups and datasets.
      \item assignment has now correct pointer target type (void *) in H5Pset_fill_value
    }
  }
}

\section{Changes in version 2.8.0}{
  \subsection{NEW FEATURES}{
    \itemize{
      \item New function h5version implemented.
      \item New low level general library functions H5open, H5close, H5garbage_collect, H5get_libversion, and H5Dset_extent implemented.
    }
  }
  \subsection{USER VISIBLE CHANGES}{
    \itemize{
      \item h5createDataset automatically uses chunking and compression.
      \item Added a warning if chunk size is equal to dimensions for large compressed datasets.
    }
  }
  \subsection{BUG FIXES}{
    \itemize{
      \item C-stack overflow when reading large fixed-length strings.
      \item error in i/o with chunksize or blocksize parameters.
      \item compiling errors due to missing int return value.
    }
  }
}

\section{Changes in version 2.6.0}{
  \subsection{NEW FEATURES}{
    \itemize{
      \item support for logical added
      \item support for reading attributes added (use read.attributes=TRUE)
      \item enabeled compression for data.frame in h5write
    }
  }
  \subsection{USER VISIBLE CHANGES}{
    \itemize{
      \item Use BiocStyles for package vignette
    }
  }
}

\section{Changes in version 2.4.0}{
  \subsection{NEW FEATURES}{
    \itemize{
      \item support for reading 64-bit integers added
      \item support for reading variable length strings added
      \item support for reading scalar objects added
    }
  }
  \subsection{USER VISIBLE CHANGES}{
    \itemize{
      \item NEWS.Rd added
      \item display of chunksize.pdf as a vignette avoided
    }
  }
}
\name{NEWS}
\title{News for Package 'rhdf5'}

\section{Changes in version 2.10.0}{
  \subsection{NEW FEATURES}{
    \itemize{
      \item Added support for HDF5 property lists.
      \item Added property list arguments to H5Dcreate and H5Dopen.
      \item New function h5readAttributes implemented that reads all HDF5 attributes of one object.
      \item New function h5version implemented.
      \item fillValue parameter added to h5createDataset.
      \item New low level general library functions H5Lcreate_external, H5Fis_hdf5, H5Fget_filesize, H5Fget_name, H5Pcreate, H5Pcopy, H5Pget_class, H5Pclose, H5Pclose_class, H5Pset_char_encoding, H5Pset_create_intermediate_group, H5Pset_chunk_cache, H5Pset_layout, H5Pset_chunk, H5Pget_chunk, H5Pset_deflate, H5Pset_fill_value, H5Pset_fill_time, H5Pset_alloc_time, H5Pequal implemented.
      \item Support for parallel Make (make -j)
    }
  }
  \subsection{USER VISIBLE CHANGES}{
    \itemize{
      \item A warning is shown in high level function (h5read, h5write and others), if an open HDF5 handle already exists for the specified filename.
    }
  }
  \subsection{BUG FIXES}{
    \itemize{
      \item Error in h5write for 0-length objects, as a consequence of automatic determining chunk size
      \item missing size parameter message in h5createDataset now correctly display
      \item checking for open file identifiers in h5read and h5ls now only searches for file names in open files, groups and datasets.
      \item assignment has now correct pointer target type (void *) in H5Pset_fill_value
    }
  }
}

\section{Changes in version 2.8.0}{
  \subsection{NEW FEATURES}{
    \itemize{
      \item New function h5version implemented.
      \item New low level general library functions H5open, H5close, H5garbage_collect, H5get_libversion, and H5Dset_extent implemented.
    }
  }
  \subsection{USER VISIBLE CHANGES}{
    \itemize{
      \item h5createDataset automatically uses chunking and compression.
      \item Added a warning if chunk size is equal to dimensions for large compressed datasets.
    }
  }
  \subsection{BUG FIXES}{
    \itemize{
      \item C-stack overflow when reading large fixed-length strings.
      \item error in i/o with chunksize or blocksize parameters.
      \item compiling errors due to missing int return value.
    }
  }
}

\section{Changes in version 2.6.0}{
  \subsection{NEW FEATURES}{
    \itemize{
      \item support for logical added
      \item support for reading attributes added (use read.attributes=TRUE)
      \item enabeled compression for data.frame in h5write
    }
  }
  \subsection{USER VISIBLE CHANGES}{
    \itemize{
      \item Use BiocStyles for package vignette
    }
  }
}

\section{Changes in version 2.4.0}{
  \subsection{NEW FEATURES}{
    \itemize{
      \item support for reading 64-bit integers added
      \item support for reading variable length strings added
      \item support for reading scalar objects added
    }
  }
  \subsection{USER VISIBLE CHANGES}{
    \itemize{
      \item NEWS.Rd added
      \item display of chunksize.pdf as a vignette avoided
    }
  }
}
\name{NEWS}
\title{News for Package 'rhdf5'}

\section{Changes in version 2.6.0}{
  \subsection{NEW FEATURES}{
    \itemize{
      \item support for logical added
      \item support for reading attributes added (use read.attributes=TRUE)
      \item enabeled compression for data.frame in h5write
    }
  }
  \subsection{USER VISIBLE CHANGES}{
    \itemize{
      \item Use BiocStyles for package vignette
    }
  }
}

\section{Changes in version 2.4.0}{
  \subsection{NEW FEATURES}{
    \itemize{
      \item support for reading 64-bit integers added
      \item support for reading variable length strings added
      \item support for reading scalar objects added
    }
  }
  \subsection{USER VISIBLE CHANGES}{
    \itemize{
      \item NEWS.Rd added
      \item display of chunksize.pdf as a vignette avoided
    }
  }
}
\name{readBeadFindMask}
\alias{readBeadFindMask}
\title{
  Read contents of a bead find mask file
}
\description{
  Reads the contents of a bead find mask file.  By default the entire bead find mask is read but there are options
  to enable loading either a rectangular sub-region or an arbitrary set of coordinates.
}
\usage{
  readBeadFindMask(
    beadFindFile,
    col=NA,
    row=NA,
    colMin=NA,
    rowMin=NA,
    colMax=NA,
    rowMax=NA
  )
}
\arguments{
  \item{beadFindFile}{
    Name of the bead find mask file to load.
  }
  \item{col,row}{
    Used to load the beadFind information for only a specific subset of wells.  col and row are a pair
    of 0-based integer vectors specifying the (col,row) coordinates of the wells to load.  If these
    options are both set to NA then the set of wells to load is determined by the values of the
    colMin,colMax,rowMin,rowMax parameters.
  }
  \item{colMin,colMax,rowMin,rowMax}{
    Used to load a rectangular sub-region of the chip.  These values are ignored if an explicit set
    of wells to load was specified by the (col,row) arguments.  All 4 values are 0-based.  If the min values
    are NA then they are replaced by 0 and if the max values are NA they are replaced by the largest possible
    coordinate.  So if all four of these options are at their default NA values the entire bead find mask
    will be read.
  }
}
\value{
  \item{beadFindMaskFile}{
    Name of the bead find mask file
  }
  \item{nCol,nRow}{
    Number of columns and rows in the associated chip
  }
  \item{col,row}{
    Integer vectors with the 0-based column and row coordinates of the wells for which bead find mask data was read.
  }
  \item{maskEmpty}{
    Boolean vector identifying the wells flagged as empty
  }
  \item{maskBead}{
    Boolean vector identifying the wells flagged as containing a bead
  }
  \item{maskLive}{
    Boolean vector identifying the wells flagged as live
  }
  \item{maskDud}{
    Boolean vector identifying the wells flagged as duds
  }
  \item{maskAmbiguous}{
    Boolean vector identifying the wells flagged as being ambiguous
  }
  \item{maskTF}{
    Boolean vector identifying the wells flagged as having a Test Fragment key signal
  }
  \item{maskLib}{
    Boolean vector identifying the wells flagged as having a library key signal
  }
  \item{maskPinned}{
    Boolean vector identifying the wells flagged as being pinned (i.e. saturated at some point in the run)
  }
  \item{maskIgnore}{
    Boolean vector identifying the wells flagged to be ignored in downstream analysis
  }
  \item{maskWashout}{
    Boolean vector identifying the wells flagged as having been washed out during the course of the run
  }
  \item{maskExclude}{
    Boolean vector identifying the wells flagged as being in the exclusion zone (not expected to yield any useful data)
  }
  \item{maskKeypass}{
    Boolean vector identifying the wells flagged as having passing key signal
  }
}
\author{
  Simon Cawley
}
\seealso{
  \code{\link{readBeadFindMaskHeader}},
  \code{\link{readWells}},
}
\name{findMixedReads}
\alias{findMixedReads}
\title{
  Estimate the (col,row) coordinates of mixed reads using results from libMixtureAnalysis()
}
\description{
  Applies the model returned from libMixutureAnalysis to all library reads to estimate the (col,row) coordinates
  of those which are mixed.
}
\usage{
  findMixedReads(
    dataDir,
    libKeySeq,
    mixResults,
    keySnrThreshold = 3
  )
}
\arguments{
  \item{dataDir}{
    The directory in which the 1.wells and bfmask.bin files are located.
  }
  \item{libKeySeq}{
    The key sequence for library fragments
  }
  \item{mixResults}{
    a list object returned from libMixtureAnalysis() which identifies which flows should be used
    for estimation of mixed reads and the mixture model to apply to them.
  }
  \item{keySnrThreshold}{
    The SNR threshold to apply to key sequences, beads with values below this will not be considered as keypass.
  }
}
\author{
Simon Cawley
}

\seealso{
  \code{\link{libMixtureAnalysis}}
}
\examples{

dataDir   <- "/data/my_analysis_directory"
libKeySeq <- "TCAG"
tfKeySeq  <- "ATCG"
tfSeq     <- "ATCGTAGCGTACATCGCGCATCTATATATCGTCAACTACGCTGAGTCGGAGACACGCAGGGATGAGATGG"

## Perform mixture analysis, use to apply to all library reads
#mixModel <- libMixtureAnalysis(dataDir, tfKeySeq, libKeySeq, tfSeq)
#mixRowCol <- findMixedReads(dataDir, libKeySeq, mixModel)
}
\name{keyStats}
\alias{keyStats}
\title{
  Compute summmary statistics based on key flows
}
\description{
  Given a matrix of raw signal intensities, keyStats() computes for each well summary metrics based on the key flows incuding the median and standard deviation of 0-mers and 1-mers in the key, as well as the total signal, noise and signal-to-noise ratio when discriminating between 1-mers and 0-mers.  
}
\usage{
  keyStats(measured,keySeq,flowOrder,sdFudge=0)
}
\arguments{
  \item{measured}{
    A matrix of signal intensities, one row per well and one column per flow.
  }
  \item{keySeq}{
    A character string representing the key sequence.
  }
  \item{flowOrder}{
    A character string representing the nucleotide flow order.
  }
  \item{sdFudge}{
    A fudge-factor to allow for avoidance of divide-by-zero issues when computing key SNR.  Any key standard deviation that is less than sdFudge is replaced by sdFudge.  Default value for sdFudge is 0.  The last 1-mer in the key is always ignored as its signal will be some mix of key and library sequence.
  }
}
\value{
  keyStats() returns a list whose elements are:
  \item{key_1_med, key_0_med}{
    Median value over the 1-mers,0-mers in the key flows.
  }
  \item{key_1_sd, key_0_sd}{
    Standard deviation over the 1-mers,0-mers in the key flows.
  }
  \item{key_sig}{
    Key signal - the difference between the median 1-mer and the median 0-mer
  }
  \item{key_sd}{
    Key standard deviation - the pooled standard deviation for the difference between 1-mers and 0-mers.  In other words, sqrt(key.0.sd^2 + key.1.sd^2)
  }
  \item{key_snr}{
    Key signal-to-noise ratio - the ratio of key_sig to key_sd.  If sdFudge is  positive then any key_sd values less than sdFudge are set to sdFudge.
  }
}
\author{
  Simon Cawley
}
\name{normalizeIonogram}
\alias{normalizeIonogram}
\title{
  Key-normalized wells data
}
\description{
  Scales all flows of each well such that the average of the key 1-mer signals is equal to 1.  This
  function assumes that all wells being normlized have the same key sequence.  It can handle being
  supplied with just a subset of all flows so long as all of the key flows are present, otherwise
  unexpected things may happen.
}
\usage{
  normalizeIonogram(
    measured,
    keySeq,
    flowOrder
  )
}
\arguments{
  \item{measured}{
    A matrix of raw signal values, a row for each well and a column for each flow.  It is fine
    to supply a subset of flows so long as all of the key flows are present, otherwise unexpected
    things may happen.
  }
  \item{keySeq}{
    The key sequence against which to normalize.
  }
  \item{flowOrder}{
    A single characer string specifying the flow order.  This string needs to be long enough to cover
    the key flows.
  }
}
\value{
  The return value is a list with a single element named normalized.
  \item{normalized}{
    A numeric vector of the key-normalized values, dimensions are the same
    as the measured vector supplied as input.
  }
}
\author{
  Simon Cawley
}
\examples{
  raw  <- matrix(rnorm(5000),nrow=100,ncol=50)
  norm <- normalizeIonogram(raw,"ATCG","TACGTACGTACG")
}
\name{flowToSeq}
\alias{flowToSeq}
\title{
  Convert flow-space values into DNA seqeunce
}
\description{
  Given a matrix of flow values and a flow order, returns the sequence under ideal circumstances - i.e. in the absence of phase errors or signal decay.
}
\usage{
  flowToSeq(flowVals,flowOrder)
  
}
\arguments{
  \item{flowVals}{
    Matrix of flow values with one row per sequence and one column per flow.  Values will be rounded to nearest integer.
  }
  \item{flowOrder}{
    A DNA string specifying the flow sequence.  Will be cycled as necessary.
  }
}
\value{
  The return value is a character vector representing the DNA sequence under ideal circumstances.
}
\author{
  Simon Cawley
}
\examples{
mySeq1     <- "TCAGCTTGTAACAGGTCAGTTACCGTCCGTCCACGCCGCCGCG"
mySeq2     <- "TCAGGCAATCAACTGGCGAAACTGGAACCGATTGTTTCGGTA"
flowOrder <- "TACG"
flowVals1  <- seqToFlow(mySeq1,flowOrder,nFlow=16)
flowVals2  <- seqToFlow(mySeq2,flowOrder,nFlow=16)
flowToSeq(rbind(flowVals1,flowVals2),flowOrder)
}
\seealso{
  \code{\link{seqToFlow}},
}
\name{readIonBam}
\alias{readIonBam}
\title{
  Read Ion Torrent BAMs with flowspace data
}
\description{
  An Ion Torrent-optimized BAM reader
}
\usage{
  read(
    bamFile,
    col=numeric(),
    row=numeric(),
    maxBases=250,
    nSample=0,
    randomSeed=1,
    readGroups=character(),
    wantMappingData=TRUE,
    maxCigarLength=100
  )
}
\arguments{
  \item{bamFile}{
    Name of the BAM file to load
  }
  \item{col,row}{
    As an alternative to returning the entire BAM file, an integer vector of 0-indexed col and row
    coordinates can be supplied to specify an arbitrary collection of wells.  If any of the wells
    that are specified are not found a warning will be issued and those that can be found are returned.
    Read coordinates are determined by parsing read names of the form "HASH:row:col".
  }
  \item{maxBases}{
    The returned bases, quality scores and flowIndex values will be truncated to this length.  The
    identity of reads that have been truncated by this limit can be established be comparing the
    return values of length and fullLength which report the trimmed and untrimmed lengths respectively.
  }
  \item{nSample}{
    As an alternative to explicitly specifying a subset of reads to load via the
    col and row arguments, a number of reads to randomly sample can be specified
    via the nSample option.
  }
  \item{randomSeed}{
    Can be used to specify the random seed when getting a random sample of
    reads.  Must be a positive integer.
  }
  \item{readGroups}{
    A character vector specifying read groups to restrict to.  To see which
    read groups are availalbe, consider using readBamReadGroups()
  }
  \item{wantgMappingData}{
    If true, alignment fields will be parsed from the BAM file
  }
  \item{maxCigarLength}{
    Sets the maximum number of elements of the cigar string that will be returned.
  }
}
\value{
  \item{nFlow}{
    The number of flows.
  }
  \item{id}{
    The read name.
  }
  \item{col,row}{
    Vectors with the 0-indexed col and row coordinates of each read.
  }
  \item{length,fullLength}{
    Vectors with number of base calls in the possibly-truncated and full-length version of each read.
  }
  \item{clipQualLeft,clipQualRight,clipAdapterLeft,clipAdapterRight}{
    Vectors with the 1-indexed quality and adapter clip positions for each read.  See SFF format documentation for more info.
  }
  \item{flowClipLeft,flowClipRight}{
    Vectors with the 0-indexed flows for the first aligned base and for the first base of any 3' adapter that was detected.
  }
  \item{flow}{
    Matrix of corrected flow signal values, one row for each read and one colum for each flow.
  }
  \item{base}{
    Vector of character strings with the called bases.  Each string is padded out with N to length maxBases as necessary.
  }
  \item{qual}{
    Matrix of per-base quality scores in Phred (-10log10) scale.  One row per read, one column per base position.
  }
  \item{flowIndex}{
    Matrix of flowIndex values mapping base calls to flows from which each base call originates.  Values are 1-indexed and
    each value is stored relative to the previous value in the array.
  }
  \item{alignFlag}{
    The alignment flag
  }
  \item{alignBase}{
    Aligned query bases (padded with dashes)
  }
  \item{alignRefID}{
    ID of the reference sequenc to which the read aligns
  }
  \item{alignPos}{
    1-based leftmost alignment position
  }
  \item{alignMapq}{
    Mapping quality
  }
  \item{alignBin}{
    Bin associated with the alignment
  }
  \item{alignCigarType}{
    Character vector, each entry is the set of cigar operation types (MIDNSHPIX=).  Vector max
    length is determined by the maxCigarLength option.
  }
  \item{alignCigarLen}{
    Integer matrix with a row for each alignment, each row describes the length of the cigar operations for the
    corresponding read.  Number of columns is limited by the maxCigarLength option
  }
  \item{header}{
    The header section of the bam file as returned by readBamHeader()
  }
}
\author{
  Simon Cawley
}
\seealso{
  \code{\link{readBeadHeader}}
}
\name{phaseFit}
\alias{phaseFit}
\title{
  Estimate phase-related parameters from flow data.
}
\description{
  Given observed flow values and the assumed-true underlying sequence, returns esimates
  of phase-related parameters for one of a few different types of model.
}
\usage{
  phaseFit(
    trueSeq,
    sig,
    flowOrder,
    cf                 = 0.0100,
    ie                 = 0.0050,
    dr                 = 0.0015,
    hpScale            = 1,
    conc               = diag(4),
    maxAdvances        = 2,
    maxIter            = 30,
    droopType          = c("ONLY_WHEN_INCORPORATING","EVERY_FLOW"),
    fitType            = c("CfIeDr","CfIeDrHpScale","HpScale","CfIeDrHpScale4","HpScale4","CfIe","NucContam","CfIe4","NucContamIe"),
    resType            = c("SQUARED","ABSOLUTE","GEMAN_MCCLURE"),
    resSummary         = c("MEAN","MEDIAN","MEAN_OF_MEDIAN"),
    ignoreHPs          = FALSE,
    flowWeight         = NULL,
    maxErr             = 1,
    extraTaps          = 0
  )
}
\arguments{
  \item{trueSeq}{
    The assumed-known sequence underlying the signal values.  Can be either a vector of DNA strings or a matrix
    of flow values, one row per read and one column per flow.
  }
  \item{sig}{
    The matrix of observed signal values, one row per read and one column per flow.
  }
  \item{flowOrder}{
    The flow cycle - for example "TACG".
  }
  \item{cf,ie,dr}{
    Estimates for cf, ie and dr.  Can be scalars, if vectors then values will be cycled over flows.
  }
  \item{hpScale}{
    HpScaling factor - incorporation signals for an HP of length h will be modeled as h*hpScale^(h-1).  Can be of length 1 or 4, in the
    case of the latter it is interpreted as a vector of per-nuc values in the order A,C,G,T.
  }
  \item{conc}{
    Estimate for the 4x4 nucleotide concentration matrix.  Column and row order is ACGT.  The value in
    row i and colum j is the amount of nucleotide j that is present when flowing nucleotide i.
    The default is to use the identity matrix.
  }
  \item{maxAdvances}{
    The maximum number of homopolymer stretches that can be extended in a single flow.
  }
  \item{maxIter}{
    The maximum number of iterations in the LevMar fit used to optimized phase parameter estimates.
  }
  \item{droopType}{
    The droop model used - can be either "ONLY_WHEN_INCORPORATING" (the default) or "EVERY_FLOW".
  }
  \item{fitType}{
    The phase model to use.  Available models are:
    \itemize{
      \item "CfIe" - fit a single carry-forward and a single incomplete extension parameter.
      \item "CfIeDr" - fit one carry-forward, one incomplete extension and one droop parameter.
      \item "CfIeDrHpScale" - fit one carry-forward, one incomplete extension, one droop and one hpScale parameter.
      \item "HpScale" - fit one hpScale parameter.
      \item "CfIeDrHpScale4" - fit one carry-forward, one incomplete extension, one droop and four nuc-specific hpScale parameters.
      \item "HpScale4" - fit four nuc-specific hpScale parameters.
      \item "NucContam" - fit the 12 off-diagonal elements of a nucleotide contamination matrix.
      \item "CfIe4" - fit a single carry-forward and 4 nucleotide-specific incomplete extension parameters.
      \item "NucContamIe" - fit a single incomplete extension parameter and the 12 off-diagonal elements of a nucleotide contamination matrix.
    }
  }
  \item{resType}{
    The type of residuals to use when evaluating the fit.  Options are:
    \itemize{
      \item "SQUARED" - r^2
      \item "ABSOLUTE" - |r|
      \item "GEMAN_MCCLURE" - x^2 / (2*(1+x^2))
    }
  }
  \item{resSummary}{
    The way transformed residuals should be summarized whene valuating the model fit.  Options are:
    \itemize{
      \item "MEAN" - the weighted mean residual across all reads and flows.
      \item "MEDIAN" - the median residual for all reads and flows with positive weight.
      \item "MEAN_OF_MEDIAN" - the average of the per-read median residuals.
    }
  }
  \item{ignoreHPs}{
    If set to true then flows involving a homopolymer stretch of more than 1 will be ignored during the fit.
  }
  \item{flowWeight}{
    A vector of weights in the range [0,1] to allow for down-weighting or ignoring certain flows during the fit.
    Length must be equal to number of flows.  Default is to set weight to 1 for every flow.
  }
  \item{maxErr}{
    The max rounding error before we start ignoring the read.  This only applies when the seuqnce is not
    set explicitly but is determined by rounding the signal.  The first time the difference between a signal
    and the nearest int is larger than this value, the weights for the rest of the read and for the 4 flows before
    are set to zero.
  }
  \item{extraTaps}{
    Controls the amount of extra flows to apply after each nuc flow.  The idea is to model situations where
    extra flows are applied to try drive to complete extension, though signal isn't actually collected on these
    flows.
  }
}
\value{
  The return value is a list with the following slots.
  \item{nIter}{
    The number of LevMar iterations to reach convergence.
  }
  \item{param.*}{
    Any fitted parameters returned will be in slots prefixed by "param."
  }
  \item{residualSummarized}{
    The optimum summarized residual value that was obtained by the LevMar fit.
  }
  \item{residualRaw}{
    A matrix of residuals (observed signal minus fitted) corresponding to the the fitted parameters.  One row
    per read and one column per flow.
  }
  \item{residualWeighted}{
    As above, but residuals are weighted as in the LevMar fit.
  }
}
\seealso{
  \code{\link{SimulateCAFIE}}, \code{\link{phaseSolve}},
}
\author{
  Simon Cawley
}
\name{readBeadFindMaskHeader}
\alias{readBeadFindMaskHeader}
\title{
  Read header section of a bead find mask file
}
\description{
  Read the header section of a bead find mask file (often named
  bfmask.bin) to determine the number of columns and rows in the
  corresponding chip.
}
\usage{
  readBeadFindMaskHeader(beadFindFile)
}
\arguments{
  \item{beadFindFile}{
    The bead find mask file to be read.
  }
}
\value{
  \item{nCol}{
    Number of columns in the associated chip
  }
  \item{nRow}{
    Number of rows in the associated chip
  }
}
\author{
  Simon Cawley
}
\seealso{
  \code{\link{readBeadFindMask}},
  \code{\link{readWells}},
}
\examples{
#bfHeader <- readBeadFindMaskHeader("/data/my_analysis_directory/bfmask.bin")
}
\name{readSFF}
\alias{readSFF}
\title{
  Read SFF (Standard Flowgram File)
}
\description{
  Reads Ion Torrent SFF files, which contain base calls and flow signal values.
}
\usage{
  readSFF(
    sffFile,
    col=numeric(),
    row=numeric(),
    maxBases=250,
    nSample=0,
    randomSeed=1
  )
}
\arguments{
  \item{sffFile}{
    Name of the SFF file to load
  }
  \item{col,row}{
    As an alternative to returning the entire SFF file, an integer vector of 0-indexed col and row
    coordinates can be supplied to specify an arbitrary collection of wells.  If any of the wells
    that are specified are not found a warning will be issued and those that can be found are returned.
  }
  \item{maxBases}{
    The returned bases, quality scores and flowIndex values will be truncated to this length.  The
    identity of reads that have been truncated by this limit can be established be comparing the
    return values of length and fullLength which report the trimmed and untrimmed lengths respectively.
  }
  \item{nSample}{
    As an alternative to explicitly specifying a subset of reads to load via the
    col and row arguments, a number of reads to randomly sample can be specified
    via the nSample option.
  }
  \item{randomSeed}{
    Can be used to specify the random seed when getting a random sample of
    reads.  Must be a positive integer.
  }
}
\value{
  \item{nFlow}{
    The number of flows.
  }
  \item{col,row}{
    Vectors with the 0-indexed col and row coordinates of each read.
  }
  \item{length,fullLength}{
    Vectors with number of base calls in the possibly-truncated and full-length version of each read.
  }
  \item{clipQualLeft,clipQualRight,clipAdapterLeft,clipAdapterRight}{
    Vectors with the 1-indexed quality and adapter clip positions for each read.  See SFF format documentation for more info.
  }
  \item{flow}{
    Matrix of corrected flow signal values, one row for each read and one colum for each flow.
  }
  \item{base}{
    Vector of character strings with the called bases.  Each string is padded out with N to length maxBases as necessary.
  }
  \item{qual}{
    Matrix of per-base quality scores in Phred (-10log10) scale.  One row per read, one column per base position.
  }
  \item{flowIndex}{
    Matrix of flowIndex values mapping base calls to flows from which each base call originates.  Values are 1-indexed and
    each value is stored relative to the previous value in the array.
  }
}
\author{
  Simon Cawley
}
\name{findDatDir}
\alias{findDatDir}
\title{
  Finds the location of the directory of dat files for a given analysis run.
}
\description{
  Given the path to a directory of analysis results, this function returns the path to the associated directory of raw DAT files.
  The mechanism is somewhat fragile - the location is parsed out of the processParameters.txt file that is expected to be located
  in the analysis directory.  For "--from-wells" analysis runs the processParameters.txt unfortunately doesn't include the information
  about the location of the raw DATs, nor in fact does any other file in the analysis directory.
}
\usage{
  findDatDir(analysisDir,paramFile="processParameters.txt",paramName="dataDirectory")
}
\arguments{
  \item{analysisDir}{
    The path to the directory of analysis results.
  }
  \item{paramFile}{
    The name of the file in the analysis directory that is expected to contain run information - "processParameters.txt" by default.
  }
  \item{paramName}{
    The prefix for the line that is expected to provide the information - "dataDirectory" by default.
  }
}
\value{
  The location of the raw results directory is returned, or NULL if it couldn't be found.
}
\author{
  Simon Cawley
}
\seealso{
  \code{\link{readDatList}},
}
\name{multreePhaser}
\alias{multreePhaser}
\title{
  Jointly estimate a DNA sequence from the flow data of multiple reads.
}
\description{
  Given observed flow values and phasing parameters for the process that generated
  them, returns an esimate of the underlying DNA sequence. Solves multiple reads
  of the same sequence jointly to obtain an estiamte of the sequence.
  Version: Sept. 17 / 2012
}
\usage{
  multreePhaser(
    signal,
    flowLimits,
    numFlows,
    flowOrders,
    PhaseParameters,
    keySeq      = "TCAG",
    basecaller  = c("treephaser-swan", "treephaser-adaptive"),
    verbose = 0
  )
}
\arguments{
  \item{signal}{
    The matrix of observed signal values, one row per read and one column per flow.
    If there are X reads per sequence, the number of rows in the matrix must be divisible 
    by X. The solver groups X consecutive reads to estimate nrow(signal)/X DNA sequences.
  }
  \item{activeUntilFlow}{
  	Vector of length nrow(signal). Contains the maximum number of flows that should
  	be used by the solver for any specific read. For a given read x, if flowLimits(x) = 0, 
  	the read will not at all be used in the solving process.
  }
  \item{numFlows}{
  	Number of flows in a given flow order. 
  	-Vector of length X, where X is the maximum number of reads per sequence.
  	-If a single value is passed, the solver assumes all X reads were obtained using the same
  	number of flows.
  }
  \item{flowOrders}{
    The flow cycles - for example "TACG".
    -Vector of length X, where X is the maximum number of reads per sequence.
    -If a single string is passed, the solver assumes all X reads were obtained using the same 
    flow cycle.
  }
  \item{PhaseParameters}{
  	Matrix of size (X,2) containing the phasing parameters. 
  	The first column contains the carry-forward values for the X reads, 
  	the second contains the incomplete extension values for the X reads.	
  }
  \item{keySeq}{
    The known key sequence at the start of the read - will be used
    for key normalization.  If specified as empty, no key normalization will
    be performed. Default: "TCAG"
  }
  \item{basecaller}{
    The variant of treePhaser to be used - options are "treephaser-swan" (the
    default) which performs sliding window normalization, and "treephaser-adaptive"
    which always restarts solving at the beginning of the read.
  }
  \item{verbose}{
  	Switch to have the solver print messages.
  	 (0) No messages (Default)
  	 (1) Messages are displayed after Solving/Normalization/Simulation steps.
  	 (10) All messages are displayed
  }
}
\value{
  The return value is a list with the following elements.
  \item{seq}{
    The estimated sequences.
  }
  \item{nBases}{
    A vector containing the number of Bases called.
  }
}
\examples{ 
	\dontrun{
   numFLows <- c(100, 50)
   flowOrders <- c("TACGTACGTCTGAGCATCGATCGATGTACAGC", "TACG")
   Phasing <- matrix(nrow=2, ncol=2)
   Phasing[1, ] <- c(0.08, 0.05)
   Phasing[2, ] <- c(0.05, 0.07)
   seq1 <- "TCAGACGGTAAGCTAGGTTAGCTTTAATCGGCGTTA"
   seq2 <- "TCAGGTATTACAGGTAGCTGATTAAAGCTCGCTAGCTAGGGATCCA"
   signal<- matrix(nrow=4, ncol=100)
   signal[1, ] <- SimulateCAFIE(seq1,flowOrders[1],Phasing[1,1],Phasing[1,2],0,numFlows[1])$sig
   signal[2, 1:50] <- SimulateCAFIE(seq1,flowOrders[2],Phasing[2,1],Phasing[2,2],0,numFlows[2])$sig
   signal[3, ] <- SimulateCAFIE(seq2,flowOrders[1],Phasing[1,1],Phasing[1,2],0,numFlows[1])$sig
   active_until <- c(100, 50, 100, 0)
   BaseCalls <- multreePhaser(signal, active_until, numFlows, flowOrders, Phasing, basecaller="treephaser-swan")
}
}
\author{
  Christian Koller
}
\name{SimulateAndSolveSeq}
\alias{SimulateAndSolveSeq}
\title{
  Simulates random base sequences given phase and noise paramters.
}
\description{
  The function SimulateAndSolveSeq Simulates random base sequences given phase and noise
  paramters. The return value is a matric of size [10, numBases+nchar(keySeq)] containing
  the fraction of reads with "row"th error occuring before or at base position "column".
}
\usage{
  SimulateAndSolveSeq <- function(
  numBases,
  numFlows,
  numWells,
  noiseSigma,
  PhaseParameters,
  flowOrder = "TACGTACGTCTGAGCATCGATCGATGTACAGC",
  keySeq = "TCAG",
  noNegativeSignal=TRUE,
  plotFigure=TRUE,
  randSeed=NA,
  diagonalStates=0
  )
}
\arguments{
  \item{numBases}{
    Number of bases in a random sequence
  }
  \item{numFlows}{
    Number of flows to be simulated
  }
  \item{numWells}{
    Number of sequences to be simulated
  }
  \item{noiseSigma}{
    Standard deviation of white gaussian noise that is added to the signal.
  }
  \item{PhaseParameters}{
    Phase parameters, vector of length 3 <CF, IE, Droop>
  }
   \item{flowOrder}{
     The flow order. Default: Samba "TACGTACGTCTGAGCATCGATCGATGTACAGC"
  }
  \item{keySeq}{
    The known key sequence at the start of the read. Default: "TCAG"
  }
  \item{noNegativeSignal}{
    If TRUE (default) no negative values can occur in the signal.
  }
  \item{plotFigure}{
    If TRUE (default) the function plots a figure.
  }
  \item{randSeed}{
    Initial seed for the random number generator.
  }
  \item{diagonalStates}{
    Switch to enable a diagonal state progression model.
  }
}
\value{
  The return value is list containing
  1) cumulativeErrorPos
    A matric of size [10, numBases+nchar(keySeq)] containing the fraction of reads 
    with "row"th error occuring before or at base position "column".
  2) meanQ17length
  3) meanQ20length
  4) meanQ30length
  5) meanQ47length
}
\examples{ 
	\dontrun{
        CumulativeErrorRate <- SimulateAndSolveSeq(250, 400, 1000, 0.08, c(0.007, 0.006, 0.001))
}
}
\author{
  Christian Koller, Nov. 12, 2012
}
\name{treePhaser}
\alias{treePhaser}
\title{
  Estimate DNA sequence from flow data.
}
\description{
  Given observed flow values and phasing parameters for the process that generated
  them, returns an esimate of the underlying DNA sequence.
}
\usage{
  treePhaser(
    signal,
    flowOrder,
    cf,
    ie,
    dr,
    keySeq      = "",
    basecaller  = c("treephaser-swan", "dp-treephaser", "treephaser-adaptive", "treephaser-solve"),
    diagonalStates=0,
    RecalModelFile="",
    RecalModelThreshold=4,
    xval=NA,
    yval=NA
  )
}
\arguments{
  \item{signal}{
    The matrix of observed signal values, one row per read and one column per flow.
  }
  \item{flowOrder}{
    The flow cycle - for example "TACG".
  }
  \item{cf,ie,dr}{
    Estimates for cf, ie and dr.
  }
  \item{keySeq}{
    The known key sequence at the start of the read - will be used
    for key normalization.  If not specified no key normalization will
    be performed.
  }
  \item{basecaller}{
    The variant of treePhaser to be used - options are 
    "treephaser-swan"      Sliding window adaptive normalization (default)
    "treephaser-adaptive"  Adaptive normalization
    "dp-treephaser"        Solve with older normalization method
    "treephaser-solve"     Solve without doing any normalization
  }
  \item{diagonalStates}{
    Switch to enable a diagonal state progression model.
  }
  \item{RecalModelFile}{
    Filename of a HP recalibration file to be loaded.
  }
  \item{RecalModelThreshold}{
    Lower (inclusive) threshold for model HP recalibration. (default 4)
  }
  \item{xval}{
    x coordinates of the wells. Required when RecalModelFile is provided.
  }
  \item{yval}{
    y coordinates of the wells. Required when RecalModelFile is provided.
  }
}
\value{
  The return value is a list with the following elements.
  \item{seq}{
    The estimated sequence.
  }
  \item{predicted}{
    The flow values predicted from the estimated sequence.
  }
  \item{residual}{
    The flow residuals (observed flow values minus predicted).
  }
  \item{hpFlow}{
    The estimated number of bases per flow.
  }
}
\seealso{
  \code{\link{SimulateCAFIE}},
}
\examples{
key   <-"TCAG"
mySeq <- paste(key,"CGCCAGGCGTTGAAGATACGCAGCGGGGCAAGCTATCCAAGGCTTCGG",sep="")
flow  <- "TACGTACGTCTGAGCATCGATCGATGTACAGC"
cf    <- 0.01
ie    <- 0.005
dr    <- 0.001
nflow <- 100
data.sim <- SimulateCAFIE(mySeq,flow,cf,ie,dr,nflow)
data.sol <- treePhaser(data.sim$sig,flow,cf,ie,dr,key)
}
\author{
  Simon Cawley
}
\name{readBamHeader}
\alias{readBamHeader}
\title{
  Read header section of a BAM file
}
\description{
  Read BAM file header
}
\usage{
  read(
    bamFile
  )
}
\arguments{
  \item{bamFile}{
    Name of the BAM file to read
  }
}
\value{
    \item{ReadGroup}{
      A list of character vectors, each vector has one entry for each read group.
      Entries are an empty string where not specified in the BAM.  Each vector
      corresponds to one of the RG tags as defined in the SAM format specification
      defined at http://samtools.sourceforge.net/SAM1.pdf
    }
    \item{Sequence}{
      A list of character vectors, each vector has one entry for each reference sequence.
      Entries are an empty string where not specified.  Each vector
      corresponds to one of the SQ tags as defined in the SAM format specification
      defined at http://samtools.sourceforge.net/SAM1.pdf
    }
}
\author{
  Simon Cawley
}
\seealso{
  \code{\link{readIonBam}}
}
\name{snrStats}
\alias{snrStats}
\title{
  Signal-to-noise analysis for wells with known sequence.
}
\description{
  Given a matrix of signal values and a true sequence associated with them, this function computes a set of summary statistics relating to typical
  signals and noise of 0-mers and 1-mers as a function of flow.  The returned data are in a form suitable for passing to other related functions
  for plotting and further analysis.

  After determining the expected ideal flow given the true sequence, the raw data are key-normalized and CAFIE-corrected.  CAFIE-estimation is
  done with findBestCafie unless the CAFIE parameters are explicitly supplied.

  The cafie-corrected signal is then analyzed separately for flows contining known 0-mers and known 1-mers.  For both 0-mers and 1-mers two robust
  linear regresisons are performed, regressing the median and sd of signal on flow, and the fitted models are part of what is returned by the
  function.
}
\usage{
snrStats(trueSeq,flowOrder,keySeq,signal,cf=NA,ie=NA,dr=NA)
}
\arguments{
  \item{trueSeq}{
    Character string representing the assumed-known sequence for all the wells being analyzed.
  }
  \item{flowOrder}{
    Character string representing the nucleotide flow order.
  }
  \item{keySeq}{
    Character string representing the key sequence for all the wells being analyzed.
  }
  \item{signal}{
    The matrix of raw signal values to be analyzed.  One row per well, one column per flow.
  }
  \item{cf,ie,dr}{
    Optionally a set of pre-computed carry-forward, incomplete extension and droop estimates can
    be supplied.  If none are supplied they are estimated on-the-fly using findBestCafie.
  }
}
\value{
  snrStats() returns a list with elements listed below.  The list is suitable for passing on to the functions snrPlotMedSig, snrPlotSdSig and oneMerSNR.
  \item{trueHP}{
    The ideal expected signal (assuming no CAFIE or droop) based on the known sequence and the flow order.
  }
  \item{medSig}{
    Vector with the median signal for each flow.
  }
  \item{sdSig}{
    Vector with the standard deviation of signal for each flow
  }
  \item{fit.med}{
    A list of two elements named "0" and "1", each of which is the robust linear model fit returned by rlm() for 0-mer and 1-mer median signals respectively.
  }
  \item{fit.sd}{
    Like fit.med, but with robust linear model fits for standard deviations of signals.
  }
  \item{cf,ie,dr}{
    The carry-forward, incomplete-extension and droop values used in CAFIE-correction.
  }
}
\author{
  Simon Cawley
}
\name{chipPlot}
\alias{chipPlot}
\title{
  Create a heatmap of data binned by (x,y) coordinates.
}
\description{
  Given a set of (x,y,z) triples, makes a heatmap of their values after doing
  some binning by (x,y) values.  Useful for making heatmaps of data by
  position on a chip, for example.
}
\usage{
  chipPlot(
    zCol,
    zRow,
    zVal,
    minCol=0,
    minRow=0,
    maxCol=NA,
    maxRow=NA,
    zlim=NA,
    nColBin=100,
    nRowBin=100,
    minBin=5,
    header="",
    histLim=NA,
    doPlot=TRUE,
    doInterpolate=FALSE,
    cex.header=1,
    summaryFunction=c("mean","median"),
    color=rgb(rep(0,256),seq(0,1,length=256),seq(1,0,length=256))
  )
}
\arguments{
  \item{zCol,zRow,zVal}{
    Vectors (of equal length) specifying the (x,y,z) triples.
    IMPORTANT: x and y values are expected to be 0-based.
  }
  \item{minCol,minRow}{
    The minimum x and y values to use when binning - zero by default.
    All x values should be greater than or equal to minCol, similar for y.
  }
  \item{maxCol,maxRow}{
    The maximum x and y values to use when binning.  If NA (the default) then
    the maximum will be set to 1 plus the max observed value.  All x values
    should be stricly less than maxCol, similar for maxRow.
  }
  \item{zlim}{
    A numeric vector of length 2 specifying lower and upper limits to
    which data will be truncated prior to plotting.  The default (if zlim
    is set to NA) is to truncate to the 2nd and 98th percentiles.
  }
  \item{nColBin,nRowBin}{
    The number of bins into which to aggregate data along X and Y axes,
    default is 100.
  }
  \item{minBin}{
    The minimum number of values per bin (default is 5).  Any bin having
    fewer than this number of values will have its entry set to NA.
  }
  \item{header}{
    Text string for plot title.  Default is no plot title.
  }
  \item{histLim}{
    A numeric vector of length 2 specifying lower and upper limits for the
    data going into the histogram plot.  If set to NA the default behaviour is
    to restrict to the range of the input image data, after applying any
    truncation implied by the zlim parameter.
  }
  \item{doPlot}{
    If set to FALSE no plot is generated, only effect is to generate
    return data.
  }
  \item{doInterpolate}{
    If set to TRUE, the binned and summarized values are bi-linearly
   interpolated back out chip-wide and returned.
  }
  \item{cex.header}{
    String expansion applied to plot title identified by header parameter.
  }
  \item{summaryFunction}{
    The function to use to summarize the data in each bin before plotting.
    Must be either "mean" (the default) or "median".
  }
  \item{color}{
    The color palette used for the image plot.  Default is a 256-step transition from blue to green.
  }
}
\value{
  \item{nColBin,nRowBin}{
    The number of column and row bins in the binned data.
  }
  \item{binnedCol,binnedRow,binnedVal}{
    The x,y and summarized aggregated z values for the binned data.
  }
  \item{sig}{
    If doInterpolate is set to TRUE, interpVal will contain z values
    after bi-linear interpolation back out to the dimensions of the input data.
  }
}
\author{
  Simon Cawley
}
\seealso{
  \code{\link{imageWithHist}}, \code{\link{formImageMatrix}}, \code{\link{bin2D}}
}
\examples{
nCol <- 200
nRow <- 100
x <- rep(0:(nCol-1),rep(nRow,nCol))
y <- rep(0:(nRow-1),nCol)
z <- rnorm(nCol*nRow,mean=x+y)
ret <- chipPlot(x,y,z,nColBin=20,nRowBin=10)
}
\name{libMixtureAnalysis}
\alias{libMixtureAnalysis}
\title{
  Analysis of library read mixtures
}
\description{
  Estimates the proportion of non-clonal wells among the library reads, using test fragment reads to estimate
  the typical locations of zero-mers and one-mers.
}
\usage{
  libMixtureAnalysis(
    dataDir,
    libKeySeq,
    plotDir = NA,
    setName = NA,
    nSample = 10000,
    densityBandwidth = 0.05,
    keySnrThreshold = 3,
    minLibNumber=5000,
    plotType = c("none", "png", "bitmap"),
    ret = TRUE
  )
}
\arguments{
  \item{dataDir}{
    The directory in which the 1.wells and bfmask.bin files are located.
  }
  \item{libKeySeq}{
    The key sequence for library fragments
  }
  \item{plotDir}{
    Directory in which plots will be created.  Optional.  If used, caller should ensure that it exists first.
  }
  \item{setName}{
    The name that will be used to refer to this dataset in plot titles.  Optional, default value is the base name of the data directory.
  }
  \item{nSample}{
    The maximum number of library and test fragment beads that will be sampled when assessing mixtures.
  }
  \item{densityBandwidth}{
    The value of the "bw" parameter that will be used in calls to density().  
  }
  \item{keySnrThreshold}{
    The SNR threshold to apply to key sequences, beads with values below this will not be considered as keypass.
  }
  \item{minLibNumber}{
    Analysis will be performed only if the number of library beads is at least this large.
  }
  \item{plotType}{
    Type of plots to produce - default value is "none".  Other options are "png" and "bitmap" which will both produce
    png files, the latter option may be requied when running in situations where there is no X running.
  }
  \item{ret}{
    If this is set to TRUE the function will return summary details and the estimated signal densities which can be used
    for testing other beads for mixing.
  }
}
\author{
  Simon Cawley
}

\seealso{
  \code{\link{findMixedReads}}
}
\examples{

dataDir   <- "/data/my_analysis_directory"
libKeySeq <- "TCAG"
tfKeySeq  <- "ATCG"
tfSeq     <- "ATCGTAGCGTACATCGCGCATCTATATATCGTCAACTACGCTGAGTCGGAGACACGCAGGGATGAGATGG"

# Perform mixture analysis and store results
#result <- libMixtureAnalysis(dataDir, tfKeySeq, libKeySeq, tfSeq)
#str(result)

# Perform mixture analysis and generate summary plots
#plotDir <- "my_plots"
#system(paste("mkdir ",plotDir))
#libMixtureAnalysis(dataDir, tfKeySeq, libKeySeq, tfSeq, plotDir=plotDir, plotType="png", ret=FALSE)
}
\name{seqToFlow}
\alias{seqToFlow}
\title{
  Convert DNA sequence into flow-space representation
}
\description{
  Given a DNA sequence and a flow order, returns the flow-space signal that would be observed under ideal circumstances - i.e. in the absence of phase errors or signal decay.
}
\usage{
  seqToFlow(sequence,flowOrder,nFlow=NA,finishAtSeqEnd=FALSE,flowOffset=0)
  
}
\arguments{
  \item{sequence}{
    The DNA string to be transformed.  Any non [ACGT] characters will be silently ignored.
  }
  \item{flowOrder}{
    A DNA string specifying the flow sequence.  If nFlow is specified the flowOrder will be cycled as necessary.
  }
  \item{nFlow}{
    The maximum number of flows to return.  If not specified, the number of flows returned will
    be equal to the length of flowOrder.
  }
  \item{finishAtSeqEnd}{
    If true, returned flow sequence will end if/when the supplied template completes,
    otherwise the returned sequence will continue out to the number of flows.
  }
  \item{flowOffset}{
    If positive, specifies a number of initial flows that will be skipped, the corresponding flow
    values will be set to NA.  Designed to work with the value stored in the ZF tag in Ion BAM files,
    returned as flowOffset by readIonBam()
  }
}
\value{
  The return value is an integer vector representing the flow values expected under ideal circumstances.
}
\author{
  Simon Cawley
}
\examples{
mySeq1     <- "TCAGCTTGTAACAGGTCAGTTACCGTCCGTCCACGCCGCCGCG"
mySeq2     <- "TCAGGCAATCAACTGGCGAAACTGGAACCGATTGTTTCGGTA"
flowOrder <- "TACG"
flowVals1  <- seqToFlow(mySeq1,flowOrder,nFlow=16)
flowVals2  <- seqToFlow(mySeq2,flowOrder,nFlow=16)
}
\seealso{
  \code{\link{flowToSeq}},
  \code{\link{readIonBam}},
}
\name{sdFromIQR}
\alias{sdFromIQR}
\title{
  Robust estimate of SD, using Inter-Quartile Range
}
\description{
  A quick-and-dirty robust estimator of standard deviation from the Inter-Quartile Range of a collection of values.
  Computes the IQR and scales it by a factor that would return an estimate of the SD if the data were Gaussian.
}
\usage{
  sdFromIQR(x,na.rm=FALSE)
}
\arguments{
  \item{x}{
    The data
  }
  \item{na.rm}{
    If set the TRUE then NA values are ignored.  Default is FALSE.
  }
}
\value{
  An estimate of the Standard Deviation
}
\author{
  Simon Cawley
}
\examples{
sdFromIQR(rnorm(1000))
}
\name{readDat}
\alias{readDat}
\title{
  Read raw Ion Torrent .dat files
}
\description{
  Reads Ion Torrent raw dat files, where one dat file corresponds to a single nucleotide flow.
  Can be used to read subsets of dat files by restricting to certain frames, wells and flows.
}
\usage{
  readDat(
    datFile,
    col=numeric(),
    row=numeric(),
    minCol=-1,
    maxCol=-1,
    minRow=-1,
    maxRow=-1,
    returnSignal=TRUE,
    returnWellMean=FALSE,
    returnWellSD=FALSE,
    returnWellLag=FALSE,
    uncompress=TRUE,
    doNormalize=FALSE,
    normStart=5,
    normEnd=20,
    XTCorrect=TRUE,
    chipType="",
    baselineMinTime=0,
    baselineMaxTime=0.7,
    loadMinTime=0,
    loadMaxTime=-1
  )
}
\arguments{
  \item{datFile}{
    Character vector with names of the dat files to load
  }
  \item{col,row}{
    As an alternative to specifying a rectangular region, an integer vector of 0-indexed col and row
    coordinates can be supplied to specify an arbitrary collection of wells.  When using this approach
    the minimum spanning rectagle of all requested wells will be read from disk.  So when the wells
    are spatially confined the read will be fast, when they are dispersed it will be slower.
  }
  \item{minCol,maxCol,minRow,maxRow}{
    Can be used to specify one or more rectangular sub-regions.  Values are 0-indexed for min and 1-indexed for max,
    so the values for max should be strictly greater than the values for min.  For multiple regions, set each to
    a vector whose length is the number of regions sought.  Setting maxCol or maxRow to -1
    leads to their being re-set to the maximum possible value.  Default is one whole-chip region.
  }
  \item{returnSignal}{
    Specifies whether or not to return a matrix of raw signal data.
  }
  \item{returnWellMean,returnWellSD,returnWellLag}{
    Specifies whether or not to return a matrix of per-well mean, SD and SD of lag-1 differences for the signal data that would be returned.
  }
  \item{uncompress}{
    Specifies whether or not to uncompress if dat is written with Variable Framerate Compression (VFC)
  }
  \item{doNormalize,normStart,normEnd}{
    Specifies if normalization should be applied - this is a subtraction from each well of the average
    value in the frame range [normStart,normEnd].  Off by default, but will be forced on if XTCorrect
    is on.
  }
  \item{XTCorrect}{
    Apply electrical cross-talk correct to undo an electrical cross-talk that occurs on the 316 and 318
    chips.  Signal between pixels in columns that are multiples of 4 apart is convolved and needs to be
    deconvolved.  Not used for 314.  If XTCorrect is true (the default) then doNormalize will be forced.
  }
  \item{chipType}{
    Explicitly set chip type, to control application of XTCorrect.  Electronic
    cross talk correction only happens if XTCorrect is TRUE and if the chip type
    is either "316" or "318".  If chipType is left as an empty string the chip
    type will be guessed from the array dimensions, but if the chip is a crop
    then the guess may be incorrect, hence the need for this option.
  }
  \item{baselineMinTime, baselineMaxTime}{
    Controls whether or not reads should be baselined.  If baselineMaxTime is greater than baselineMinTime
    then each well will have the average of the flows completely contained within the timeframe subtracted.
    The average is weighted by the duration spanned by each frame used.
  }
  \item{loadMinTime,loadMaxTime}{
    Specify times, in seconds, to determine which frames are returned.
    If loadMaxTime is less than loadMinTime then all frames are returned.
  }
}
\value{
  \item{datFile}{
    Character vector with the names of the dat files loaded.
  }
  \item{nCol,nRow,nFrame,nFlow}{
    The number of columns, rows, frames and flows represented in the full dat file.
    If a subset of wells or frames has been requested nCol, nRow and nFrame will still
    describe the full dat file.  nFrame will describe the number of compressed frames
    in a dat file written with Variable Framerate Compression (VFC), unless uncompress
    is set to true in which case nFrame will describe the number of uncompressed frames.
  }
  \item{col,row}{
    Integer vectors with the 0-indexed coordinates of the wells returned.
  }
  \item{frameStart,frameEnd}{
    Vector specifying the frame start/end time in seconds.  If multiple dat files are read this vector repeats such that its length is equal to the total number of frames loaded across all flows.
  }
  \item{signal}{
    Numeric matrix of the raw signal data.  This is only returned if returnSignal is TRUE.  One column per frame loaded and one row per well loaded.  In the case of multiple dats the number of columns is equal to the number of dats times the number of frames loaded per dat.
  }
  \item{wellMean,wellSD}{
    Numeric matrices of the mean and SD of the signal data that would be returned.  Only returned if returnWellMean
    and returnWellSD are true, respectively.  One row per well and one column per dat.
  }
}
\author{
  Simon Cawley
}
\seealso{
  \code{\link{readDatCollection}},
}
\name{readDatCollection}
\alias{readDatCollection}
\title{
  Read groups of raw Ion Torrent .dat files
}
\description{
  Reads groups of Ion Torrent raw dat files, where one dat file corresponds to a single nucleotide flow.
  Can be used to read subsets of data by restricting to certain frames or wells.
}
\usage{
  readDatCollection(
    datDir=NA,
    analysisDir=NA,
    minFlow=1,
    maxFlow=-1,
    col=numeric(),
    row=numeric(),
    minCol=0,
    maxCol=-1,
    minRow=0,
    maxRow=-1,
    returnSignal=TRUE,
    returnWellMean=FALSE,
    returnWellSD=FALSE,
    uncompress=TRUE,
    baselineMinTime=0,
    baselineMaxTime=0.7,
    loadMinTime=0,
    loadMaxTime=-1
  )
}
\arguments{
  \item{datDir,analysisDir}{
    Exactly one of these must be specified.  If analysisDir is set then datDir, the location of the dat files, will be determined automatically.
    Otherwise the dat files are assumed to be located at datDir.
  }
  \item{minFlow,maxFlow}{
    Specify the first and last flows that should be loaded (1-indexed).  Setting maxFlow to -1 leads
    to it being re-set to the maximum available flow.
  }
  \item{col,row}{
    As an alternative to specifying a rectangular region, an integer vector of 0-indexed col and row
    coordinates can be supplied to specify an arbitrary collection of wells.
  }
  \item{minCol,maxCol,minRow,maxRow}{
    Can be used to specify a rectangular sub-region.  Values are 0-indexed.  Setting maxCol or maxRow to -1
    leads to their being re-set to the maximum possible value.
  }
  \item{baselineMinTime, baselineMaxTime}{
    Specify the time interval to be used for determination of well-specific baseline to subtract.
  }
  \item{returnSignal, returnWellMean, returnWellSD}{
    Specify what kind of data are returned.  If returnSignal is true the per-frame data are returned for
    all frames in the requested time interval.  If returnWellMean and returnWellSD the per-well mean and
    standard deviation of signal is returned.
  }
  \item{uncompress}{
    If true and if the data are stored with Variable Framerate Compression, the data are interpolated back out to
    all frames, at the cost of speed and memory.  If false, compressed data are returned as stored, which can be
    significantly faster and more compact.
  }
  \item{loadMinTime,loadMaxTime}{
    Specify the time interval for which to return data.  Setting loadMaxTime to -1 leads
    to it being re-set to the maximum possible value
  }
}
\value{
  The returned value is a list with the same entries as those returned by readDat (which this function wraps).
  The following additional entries are also set:
  \item{flow}{
    A vector specifying the flows that were loaded.
  }
}
\author{
  Simon Cawley
}
\seealso{
  \code{\link{readDat}},
}
\name{formImageMatrix}
\alias{formImageMatrix}
\title{
  Transform (x,y,z) triples into a matrix suitable for supplying to image() or imageWithHist()
}
\description{
  Given a set of (x,y,z) triples along with the size of the x and y
  dimensions to which the data belong, returns a matrix of z values
  with number of columns equal to the x dimension and number of rows equal
  to the y dimension.

  Matrix entries are NA where there is no corresponding (x,y) value.  The
  matrix can be supplied to \code{\link{image}} or \code{\link{imageWithHist}}
  for plotting.
}
\usage{
  formImageMatrix(
    x,
    y,
    z,
    maxX,
    maxY
  )
}
\arguments{
  \item{x,y,z}{
    Vectors (of equal length) specifying the (x,y,z) triples.

    IMPORTANT: x and y values are expected to be 0-based.
  }
  \item{maxX,maxY}{
    The number of columns (maxX) and rows (maxY) in the returned matrix.
    All x values are expected to be strictly less than maxX and all y
    values strictly less than maxY.
  }
}
\author{
  Simon Cawley
}
\seealso{
  \code{\link{imageWithHist}}, \code{\link{bin2D}}, \code{\link{chipPlot}}
}
\examples{
nRow <- 20
nCol <- 50
x <- rep(0:(nCol-1),rep(nRow,nCol))
y <- rep(0:(nRow-1),nCol)
z <- formImageMatrix(x,y,rnorm(nCol*nRow,mean=x+y),maxX=nCol,maxY=nRow)
imageWithHist(z,header="Example",zlim=range(z))
}
\name{readDatList}
\alias{readDatList}
\title{
  Returns the list of all DAT files and associated flow numbers for a run.
}
\description{
  Given a directory which may contain raw DAT files for a run, returns the names of any DAT files and associated flow numbers.
  Searches for DAT files by finding files whose name matches the expected regular expression acq_####.dat
}
\usage{
  readDatList(datDir)
}
\arguments{
  \item{datDir}{
    The directory to search.
  }
}
\value{
  The return value is a list with two entries
  \item{datFiles}{
    A character vector specifying the dat file names.
  }
  \item{datFlows}{
    A numeric vector specifying the flow number for each dat file.  Numbering starts at 1, so acq_0000.dat would be flow 1, and so on.
  }
}
\author{
  Simon Cawley
}
\seealso{
  \code{\link{findDatDir}},
}
\name{bkgModel}
\alias{bkgModel}
\title{
  Fit the background model to raw dat data
}
\description{
  Fits the background model to raw Ion Torrent dat files to produce estimates of
  incorporation signals, along with a lot of additional parameters and diagnostics.
}
\usage{
  bkgModel(
    dat,
    bkgWell,
    fitWell,
    maxFlow=NA,
    sigma_guess=2.5,
    t0_guess=23,
    tau_guess=25.0,
    dntp_uM=50.0
  )
}
\arguments{
  \item{dat}{
    A list with dat data, such as can be obtained with \code{\link{readDatCollection}}.  The list must contain the following names elements
    \itemize{
      \item signal: A matrix of signal values, each row is a well and each column is a frame.  Flows are stored consecutively, so the number of columns should be equal to the number of frames times the number of flows.
      \item nFlow: The total number of flows.
      \item nFrame: The number of frames per flow.
    }
  }
  \item{bkgWell}{
    Used to specify the rows corresponding to the empty wells in the dat$signal matrix.
    Can be either a vector of Booleans of length equal to the number of wells, or a vector
    of integers in the range [1,nrow(dat$signal)].
  }
  \item{fitWell}{
    Used to specify the wells to which the background model should be applied.  As with bkgWell,
    this vector can be Boolean or integer.
  }
  \item{maxFlow}{
    Limit the fit to flows less than or equal to maxFlow (1-indexed).  By default all available
    flows will be fit.
  }
  \item{sigma_guess}{
    Initial sigma_guess estimate supplied to the background model.
  }
  \item{t0_guess}{
    Initial t0_guess estimate supplied to the background model.
  }
  \item{tau_guess}{
    Initial tau_guess estimate supplied to the background model.
  }
  \item{dntp_uM}{
    Initial dntp_uM estimate supplied to the background model.
  }
}
\value{
  The return value is a list with estimates and parameters from the model fit.  Full documentation
  of all the values is beyond the current scope of this man page, but here's a start:
  \item{sig}{
    Matrix of estimated incorporation signals, dimension is (nWell,nFlow).
  }
  \item{bkg}{
    Vector with estimated background trace supplied to the background model, length is (nFrame*nFlow).
  }
  \item{nFitFrame}{
    The number of frames modeled by the background model.  The model compresses some frames to the
    fitted background and signal that it returns have fewer frames than are passed in.
  }
  \item{fitFg}{
    Matrix of fitted dat-level data, dimension is (nWell,nFitFrame).
  }
  \item{fitBg}{
    Vector with fitted background trace returned by the background model, length is (nFitFrame*nFlow)
  }
  \item{offset}{
    Matrix of fitted offset values, dimension is (nWell*nFlow)
  }
  \item{tau}{
    Matrix of fitted tau values, dimension is (nWell*4) - i.e. one per nucleotide per well
  }
  \item{nFlowBatch}{
    The model is fit in batches, typically 8 flows per batch.  nFlowBatch is the number of flow batches that were fitted.
  }
  \item{R,dt,P,gain,tshift}{
    Matrices of fitted parameters, dimensions are (nWell,nFlowBatch) - i.e. one per well per flow batch.
  }
  \item{krate,d,mR,oR}{
    Matrices of fitted parameters, dimensions are (4,nFlowBatch) - i.e. one per nucleotide per flow batch.
  }
  \item{sigma}{
    Vector of fitted sigma values, length is nFlowBatch.
  }
  \item{sens,C}{
    Fitted scalar parameter values.
  }
}
\author{
  Simon Cawley
}
\seealso{
  \code{\link{readDatCollection}},
}
\examples{\dontrun{
analysisDir <- "/results2/analysis/output/IonWest/B6_60_trunk_r5928_SEC_2118"
minCol <- 300
maxCol <- 349
minRow <- 300
maxRow <- 349
maxFlow <- 16

# Read 16 flows from raw dats
dat <- readDatCollection(analysisDir=analysisDir,minCol=minCol,maxCol=maxCol,minRow=minRow,maxRow=maxRow,maxFlow=maxFlow)

# Read the beadfind mask
bfFile   <- sprintf("\%s/bfmask.bin",analysisDir)
bf   <- readBeadFindMask(bfFile)

# Determine the empty wells in the slice that was loaded
regionMask <- bf$col >= minCol & bf$col <= maxCol & bf$row >= minRow & bf$row <= maxRow
emptyMask <- bf$maskBead[regionMask]==0

# Fit the background model on the first 20 non-empty wells
bgFit <- bkgModel(dat,emptyMask,which(!emptyMask)[1:20])

# Insepct what's avaialble to play with in the returned object
str(bgFit)
}}
\name{bin2D}
\alias{bin2D}
\title{
  Aggregage (x,y,z) triples into bins on (x,y) dimensions.  Can be used
  to bin z data by (x,y) values for plotting as a heatmap at coarser
  resolution than the original data.
}
\description{
  Given (x,y,z) triples, returns a new set of (x,y,z) triples binned by (x,y)
}
\usage{
  bin2D(
    x,
    y,
    z,
    minX=0,
    minY=0,
    maxX=NA,
    maxY=NA,
    nBinX=100,
    nBinY=100,
    minBin=1
  )
}
\arguments{
  \item{x,y,z}{
    Vectors (of equal length) specifying the (x,y,z) triples.
  }
  \item{minX,minY}{
    The minimum x and y values to use when binning - zero by default.
    All x values should be greater than or equal to minX, similar for y.
  }
  \item{maxX,maxY}{
    The maximum x and y values to use when binning.  If NA (the default) then
    the maximum will be set to 1 plus the max observed value.  All x values
    should be stricly less than maxX, similar for maxY.
  }
  \item{nBinX,nBinY}{
    The number of (equal-width) bins into which to split the x and y data.
  }
  \item{minBin}{
    The minum number of observations required for each bin. If the actual
    number of data is less than this minum the bin in question is set to NA.
  }
}
\author{
  Simon Cawley
}
\seealso{
  \code{\link{imageWithHist}}, \code{\link{formImageMatrix}}, \code{\link{chipPlot}}
}
\examples{
# Example of a heatmap
nRow <- 20
nCol <- 50
x <- rep(0:(nCol-1),rep(nRow,nCol))
y <- rep(0:(nRow-1),nCol)
z <- rnorm(nCol*nRow,mean=x+y)
zMat <- formImageMatrix(x,y,z,maxX=nCol,maxY=nRow)
imageWithHist(zMat,header="Fine")

# Make it coarser
nRowCoarse <- nRow/5
nColCoarse <- nCol/5
coarse <- bin2D(x,y,z,nBinX=nColCoarse,nBinY=nRowCoarse)
zSummary <- unlist(lapply(z,mean))
zMatCoarse <- formImageMatrix(coarse$x,coarse$y,zSummary,maxX=nColCoarse,maxY=nRowCoarse)
imageWithHist(zMatCoarse,header="Coarse")
}
\name{phaseSolve}
\alias{phaseSolve}
\title{
  Estimate DNA sequence from flow data.
}
\description{
  Given observed flow values and phasing parameters for the process that generated
  them, returns an esimate of the underlying DNA sequence.
}
\usage{
  phaseSolve(
    signal,
    flowOrder,
    cf,
    ie,
    dr,
    hpScale              = 1,
    conc                 = diag(4),
    droopType            = c("ONLY_WHEN_INCORPORATING","EVERY_FLOW"),
    maxAdvances          = 2,
    nIterations          = 3,
    residualScale        = TRUE,
    residualScaleMinFlow = -1,
    residualScaleMaxFlow = -1,
    extraTaps            = 0
  )
}
\arguments{
  \item{signal}{
    The matrix of observed signal values, one row per read and one column per flow.
  }
  \item{flowOrder}{
    The flow cycle - for example "TACG".
  }
  \item{cf,ie,dr}{
    Estimates for cf, ie and dr.  Can be scalars, if vectors then values will be cycled over flows.
  }
  \item{hpScale}{
    HpScaling factor - incorporation signals for an HP of length h will be modeled as h*hpScale^(h-1).  Can be of length 1 or 4, in the
    case of the latter it is interpreted as a vector of per-nuc values in the order A,C,G,T.
  }
  \item{conc}{
    The 4x4 nucleotide concentration matrix.  Column and row order is ACGT.  The value in
    row i and colum j is the amount of nucleotide j that is present when flowing nucleotide i.
    The default is to use the identity matrix.
  }
  \item{droopType}{
    The droop model used - can be either "ONLY_WHEN_INCORPORATING" (the default) or "EVERY_FLOW".
  }
  \item{maxAdvances}{
    The maximum number of homopolymer stretches that can be extended in a single flow.
  }
  \item{nIterations}{
    The maximum number of iterations of seqeunce estimation.
  }
  \item{residualScale}{
    if true, then enables signal rescaling based on residuals after each iteration of sequence estimation.
  }
  \item{residualScaleMinFlow,residualScaleMaxFlow}{
    The first and last flows to use for residual scaling.  0-based.  If set to -1 (the default) then
    whatever is the current default in the underlying C++ code will be used.
  }
  \item{extraTaps}{
    Controls the amount of extra flows to apply after each nuc flow.  The idea is to model situations where
    extra flows are applied to try drive to complete extension, though signal isn't actually collected on these
    flows.
  }
}
\value{
  The return value is a list with the following elements.
  \item{seq}{
    The estimated sequence.
  }
  \item{predicted}{
    The flow values predicted from the estimated sequence.
  }
  \item{residual}{
    The flow residuals (observed flow values minus predicted).
  }
  \item{hpFlow}{
    The estimated number of bases per flow.
  }
  \item{multiplier}{
    If residualScale is set to TRUE then the multipliers used for each iteration are returned in this slot.  The
    multiplier used are returned in reverse order - the multiplier used for the last iteration is returned first.
  }
}
\seealso{
  \code{\link{SimulateCAFIE}}, \code{\link{phaseFit}},
}
\author{
  Simon Cawley
}
\name{tnormalizeRead}
\alias{normalizeRead}
\title{
  Normalize reads.
}
\description{
  Normalize (key normalized) reads w.r.t a predicted signal.
}
\usage{
normalizeRead(
		signal,
		prediction,
		method=c("adaptive", "gain", "pid"),
		windowSize=0,
		numSteps=0,
		startFlow=0,
        endFlow=0
)
}
\arguments{
  \item{signal}{
    The matrix of observed signal values, one row per read and one column per flow.
  }
  \item{prediction}{
    The matrix of predicted signal values, one row per read and one column per flow.
  }
  \item{method}{
    Selecting the normalization method to be used:
    adaptive: Fitting a windowed additive and multiplicative offset.
    gain:     Fitting a multiplicative offset between startFlow and EndFlow.
    pid:      Fitting a multiplicative offset through a PID controller.
  }
  \item{windowSize}{
    Window size for adaptive normalization. If zero, the built-in default of the c++ code is used.
  }
  \item{numSteps}{
    The number of steps for adaptive normalization. One step computes one window 
    and interpolates between the last and current window. If 0, the number of
    steps is set to num_flows/windowSize.
  }
  \item{startFlow}{
    Flow at which to start gain of PID normalization.
  }
  \item{endFlow}{
    Flow at which to end gain or PID normalization. (Open interval)
  }
}
\value{
  The return value is a list with the following elements.
  \item{method}{
    String indicating the normalization method used.
  }
  \item{normalized}{
    Matrix containing the normalized signal values.
  }
}
\author{
  Christian Koller
}
\name{readWells}
\alias{readWells}
\title{
  Read signal data from a .wells file
}
\description{
  Reads signal data from a .wells file, offering various options to load specific subsets of information
  so that a manageable-sized chunk can be loaded into memory.  By default this function will try to read
  all flows of all wells, however in many cases that would involve reading unmanageable amounts of data
  into memory.  As a result it is likely that the most common use of this function will involve one or
  more of the flow,col,row,colMin,rowMin,colMax,rowMax options.
}
\usage{
  readWells(
    wellPath,
    col=NA,
    row=NA,
    bfMaskFile=NA,
    bfMaskFileName="bfmask.bin",
    colMin=NA,
    rowMin=NA,
    colMax=NA,
    rowMax=NA,
    ignoreBfMaskFile=FALSE,
    nCol=NA,
    nRow=NA,
    flow=numeric()
  )
}
\arguments{
  \item{wellPath}{
    Wells File to read.
  } 
  \item{flow}{
    Integer 1-based vector to enable loading signal data for only a specific subset of flows.  If this
    vector has zero-length (the default) then all available flows will be loaded.
  }
  \item{col,row}{
    Used to load signal data for only a specific subset of wells.  col and row are a pair
    of 0-based integer vectors specifying the (col,row) coordinates of the wells to load.  If these
    options are both set to NA then the set of wells to load is determined by the values of the
    colMin,colMax,rowMin,rowMax parameters.
  }
  \item{bfMaskFile}{
    The name of a bead find mask file from which to load up bead find-related information for the wells
    being read.  If this value is NA (the default) and if ignoreBfMaskFile is TRUE (the default) then
    the directory in which the wells file is located will be searched for a file named bfMaskFileName
    to load for beadFindMask information.
  }
  \item{bfMaskFileName}{
    The name of the bead find mask file to look for in the directory containing the wells file, for
    cases where the caller isn't explicitly specifying the bead find mask filename.
  }
  \item{colMin,colMax,rowMin,rowMax}{
    Used to load a rectangular sub-region of the chip.  These values are ignored if an explicit set
    of wells to load was specified by the (col,row) arguments.  All 4 values are 0-based.  If the min values
    are NA then they are replaced by 0 and if the max values are NA they are replaced by the largest possible
    coordinate.  So if all four of these options are at their default NA values the entire bead find mask
    will be read.
  }
  \item{ignoreBfMaskFile}{
    By default this is set to FALSE, but if it is TRUE then there will be no attempt to load up
    bead find mask information.  In that event, the return value will contain no bead find mask
    information.  Note that the bead find mask file is where the information on the number of 
    rows and columns is obtained (since they are not present in the wells file), so if ignoreBfMaskFile
    is set to true then nCol and nRow must be set.
  }
  \item{nCol,nRow}{
    If no bead find information is loaded then the number of columns and rows in the chip must be
    explicitly provided in this option.
  }
}
\value{
  The return value is a list with elements described below:
  \item{wellFile}{
    The name of the wells file from which data is loaded.
  }
  \item{nCol,nRow,nFlow}{
    The total number of columns, rows and flows represented in the wells file.
  }
  \item{nLoaded}{
    The total number of wells that have been read into memory.
  }
  \item{flow}{
    An integer vector (1-based) specifing the flows that were read into memory.
  }
  \item{flowBase}{
    A character vector specifying the nucleotide flowed in each one of the flows that were loaded.
  }
  \item{flowOrder}{
    A character string specifying the nucleotides flowed in the complete set of flows represented in the wells file
    (regardless of which flows have actually been loaded into memory)
  }
  \item{col,row}{
    Integer vectors (0-based) specifying the column and row coordinates of the wells that were loaded into memory.
  }
  \item{rank}{
    A numeric vector containing the rank values from each of the wells that were loaded into memory.
  }
  \item{beadFindMaskFile}{
    The bead find mask file that was read when reading the wells file (empty if none was read).
  }
  \item{mask}{
    A list of boolean vectors describing the wells that were loaded into memory.  If no bead find mask
    file was read th list will be empty.  Otherwise the list will contain entries as named below, for
    detail on what each represents see \code{\link{readBeadFindMask}}.
      empty
      bead
      live
      dud
      ambiguous
      tf
      lib
      pinned
      ignore
      washout
  }
  \item{signal}{
    A matrix with the raw signal data from the wells file.  There is a column for each flow and a row for each well.
  }
}
\author{
  Simon Cawley
}
\seealso{
  \code{\link{readBeadFindMask}},
  \code{\link{readBeadFindMaskHeader}},
}
\examples{
# x1 <- readWells("/mydata/1.wells",flows=1:10)
# x2 <- readWells("/mydata/1.wells",col=0:4,row=0:4)
}
\name{calcHypothesesDistances}
\alias{calcHypothesesDistances}
\title{
  Calculate the squared distances betweeen hypotheses base sequences.
}
\description{
  The function calculates the squared distances betweeen the predicted ionograms of 
  hypothesis base sequences to normalized observed sequences. Bases before startFlow and 
  after the end of the hypothesis are filled in by the solver based on the signal.
  Only flows are taken into where the predictions of hypotheses i, i in {1,...,N},
  differ from the predictions of hypothesis 0 by 0.05 or more.
}
\usage{
 calcDistances <- function(
  signal,
  cf,
  ie,
  dr,
  flowOrder,
  hypotheses,
  startFlow = 0,
  normalize = 0,
  verbose = 0
)
}
\arguments{
  \item{signal}{
    The vector of (key) normalized measured signal values.
    Adaptive normalization is performed for every hypothesis sequence and the
    resulting normalized values are returned in field "Normalized"
  }
  \item{flowOrder}{
    The flow cycle - for example "TACG".
  }
  \item{cf,ie,dr}{
    Estimates for cf, ie and dr.
  }
  \item{hypotheses}{
    String Vector containing the base space hypotheses.
  }
  \item{startFlow}{
    Flow which corresponds to the start of the hypotheses. Flows up to startflow
    are solved and a corresponding base prefix is attached to the hypotheses.
  }
  \item{normalize}{
  	If normalize>0, the signal values are adaptively normalized w.r.t the hypothesis prediction.ß
  }
  \item{verbose}{
  	If verbose>0, the functions prints out some text.
  }
}
\value{
  The return value is a list with the following elements.
  \item{DistanceObserved}{
    Squared distance of the normalized signal and the prediction for the hypotheses.
  }
  \item{DistanceHypotheses}{
    Squared distance of the prediction for hypotheses zero and i.
  }
  \item{Predictions}{
    Matrix of prediced signal values given the base seqences.
    Rows: different hypotheses Columns: Predicted signal for flows
  }
  \item{Normalized}{
  	Normalized measurements.
  	- Equal to input signal if normalize==0
  	- Adaptively normalized w.r.t hypothesis prediction if normaize>0
  	Rows: different hypotheses Columns: Normalized flow signals
  }  
  
}
\seealso{
  torrentR package
}
\examples{
\dontrun{
	key   <-"TCAG"
	mySeq <- paste(key, "GGCGCCAGGCGTTGAAGATACGCAGCGGGGCAAGCTATCCCCAAGGCTTCGG", sep="")
	Hyp1 <- "GGCGCCAGGCGTTGAAGATACGCAGCGGGGCAAGCTATCC"
	Hyp2 <- "GGCGCCAGGCGTTGAAGATACGCAGCGGGCAAGCTATCC"
	Hyp3 <- "GGCGCCAGGCGTTGAAGATACGCAGCGGGGGCAAGCTATCC"
	flow  <- "TACGTACGTCTGAGCATCGATCGATGTACAGC"
	startFlow <- 7
	cf    <- 0.01
	ie    <- 0.005
	dr    <- 0.001
	nflow <- 100
	signal <- SimulateCAFIE(mySeq,flow,cf,ie,dr,nflow)
	vals <- calcHypothesesDistances(signal$sig,cf,ie,dr,flow,c(Hyp1, Hyp2, Hyp3), startFlow,1)
}
}
\author{
  CK
}
\name{imageWithHist}
\alias{imageWithHist}
\title{
  Create a heatmap with an associated histogram.
}
\description{
  A wrapper around the image() function which additionally creates a histogram
  giving a sense of the distribution of the values forming the heatmap.
}
\usage{
  imageWithHist (
    z,
    zlim=NA,
    header="",
    histLim=NA,
    xaxt="n",
    yaxt="n",
    nHistBar=100,
    cex.header=1,
    col=rgb(rep(0,256),seq(0,1,length=256),seq(1,0,length=256)),
    ...
  )
}
\arguments{
  \item{z}{
    The data to be passed along to the image() function, modulo whatever
    truncation is applied due to the setting of the zlim parameter.
  }
  \item{zlim}{
    A numeric vector of length 2 specifying lower and upper limits to
    which data will be truncated prior to plotting.  The default (if zlim
    is set to NA) is to truncate to the 2nd and 98th percentiles.
  }
  \item{header}{
    Text string for plot title.  Default is no plot title.
  }
  \item{histLim}{
    A numeric vector of length 2 specifying lower and upper limits for the
    data going into the histogram plot.  If set to NA the default behaviour is
    to restrict to the range of the input image data, after applying any
    truncation implied by the zlim parameter.
  }
  \item{xaxt,yaxt}{
    The axis type for the x and y axes of the image plot.  See par() for
    more info.
  }
  \item{nHistBar}{
    The number of bars in the histogram - default is 100.
  }
  \item{cex.header}{
    String expansion applied to plot title identified by header parameter.
  }
  \item{col}{
    The color palette used for the image plot.  Default is a 256-step transition from blue to green.
  }
  \item{...}{
    Any additional options get passed through to the call to image()
  }
}
\author{
  Simon Cawley
}
\seealso{
  \code{\link{formImageMatrix}}, \code{\link{bin2D}}, \code{\link{chipPlot}}
}
\examples{
nRow <- 20
nCol <- 50
x <- rep(0:(nCol-1),rep(nRow,nCol))
y <- rep(0:(nRow-1),nCol)
z <- formImageMatrix(x,y,rnorm(nCol*nRow,mean=x+y),maxX=nCol,maxY=nRow)
imageWithHist(z,header="Example",zlim=range(z))
}
\name{torrentR_version}
\alias{torrentR_version}
\title{
  Returns the version of the torrentR package
}
\description{
  Returns the version of the torrentR package
}
\usage{
  torrentR_version()
}
\value{
  Returns a character string specifying the torrentR version
}
\author{
  Simon Cawley
}
\seealso{
  \code{\link{plot_version}},
}
\examples{
  x <- torrentR_version()
}
\name{writeFASTQ}
\alias{writeFASTQ}
\title{
  Writes a vetor of sequences into a .fastq file
}
\description{
  Writes a string vector of sequences into a .fastq file}
\usage{
  writeFASTQ <- function (
    fileNamePath,
    sequences,
    qualityValues = NA,
    wellRow,
    wellColumn,
    keySeq="TCAG",
    keyPassFilter=TRUE,
    appendFile=FALSE
)
}
\arguments{
  \item{fileNamePath}{
    Name base and path where the .fastq file should be saved
  }
  \item{sequences}{
  	Vector of DNA sequences
  }
  \item{qualityValues}{
    Vector of quality value strings (same length as sequences)
    default: NA which creates dummy quality strings in the fastq file
  }
  \item{wellRow}{
    Vector of row indices of wells that created the sequences.
  }
  \item{wellColumn}{
  	Vector of column indices of wells that created the sequences.
  }
  \item{keySeq}{
    The known key sequence at the start of the read. Default: "TCAG"
  }
  \item{keyPassFilter}{
    If TRUE (default) sequences whose first called bases do not correspond to the
    keySeq are not written to the fastq file.
  }
  \item{appendFile}{
  	If TRUE, the function does not overwrite an existing file but appends the
  	information to the end of the file. Default: FALSE
  }
}
\value{
  The return value is a vector with three elements.
  1) number of sequences at the input
  2) number of sequences that had the correct key sequence at the beginning
  3) number of sequences that failed key-pass or were too short (length 0 after removal of key)
}
\examples{ 
	\dontrun{		
	Seq = c("TCAGACGGTAAGCTAGGTTAGCTTTAATCGGCGTTA", "TCAGGTATTACAGGTAGCTGATTAAAGCTCGCTAGCTAGGGATCCA")
	logVec <- writeFASTQ("MyFastq", Seq, NA, c(0,1), c(1,3))
}
}
\author{
  Christian Koller
}
\name{plot_version}
\alias{plot_version}
\title{
  Add the torrentR version number to a plot
}
\description{
  Makes a call to mtext to write the torrentR version id to the current
  plotting device.
}
\usage{
  plot_version(
    side=1,
    line=4,
    adj=1,
    ...
  )
}
\arguments{
  \item{side}{
    The side of the plot to which the version string is added (1=bottom, 2=left, 3=top, 4=right).
  }
  \item{line}{
    The margin line to which text is added.  Starts at 0, counting outward.
  }
  \item{adj}{
    Text justification.  For strings parallel to the axes, 'adj = 0' means
    left or bottom alignment, and 'adj = 1' means right or top alignment.
  }
  \item{...}{
    Any additional arguments are passed right through to mtext().
  }
}
\author{
  Simon Cawley
}
\seealso{
  \code{\link{torrentR_version}},
}
\examples{
  plot(rnorm(100))
  plot_version()
}
\name{RepeatSeqBasecalling}
\alias{RepeatSeqBasecalling}
\title{
  Generates Basecalls from multiple wells files corresponding to the same chip, taking advantage of
  multiple wells file jointly.
}
\description{
  The function generates Basecalls from multiple wells files corresponding to the same chip. Uses 3
  input files per run, 1.wells, analysis.bfmask.bin, and BaseCaller.json. 
  Assume x in {1,2,3}, then naming of the files should be
  1.x.wells
  analysis.bfmask.x.bin
  BaseCaller.x.json
  }
\usage{
  RepeatSeqBasecalling <- function(
  DataFolder,
  numRuns,
  chipRegion=NA,
  unionOfReads=FALSE
  combinations=NA
  outputFolder="."
)
}
\arguments{
  \item{DataFolder}{
    Path where the input files are stored (without last slash)
  }
  \item{numRuns}{
  	Number of different runs in DataFolder. (Range of parameter x in Description)
  }
  \item{chipRegions}{
    A chip is divided in 8x8 regions, numbered 0-63. chipRegions is a vector input which
    specifies which regions should be processed. Default: NA -> all regions are processed
  }
  \item{combinations}{
    Logical vector / matrix with rows of length numRuns, indicating which runs should
    be called jointly. Default: All runs are called jointly.
    Please put the case that inolves the most runs in row 1. 
  }
  \item{unionOfReads}{
  	Logical switch that determines whether the union or intersection set of reads
  	across the runs is called.
  }
  \item{outputFolder}{
  	Folder where output files are written. Default: current directory
  }
}
\value{
  There is no return value, but 2 files are generated per row in combinations
  1) BaseCalls.x.y.z.fastq, containing called sequences in fastq format
  2) log.txt, number of called sequences, sequences that key-pass, and fail per region
}
\examples{ 
	\dontrun{
	combination <- matrix(nrow=3, ncol=3)
	combination[1, ] <- c(TRUE, TRUE, TRUE)
	combination[2, ] <- c(TRUE, FALSE, TRUE)
	combination[3, ] <- c(TRUE, FALSE, FALSE)
	
	RepeatSeqBasecalling("./", 3, c(18, 45), combinations=combination)
    }
}
\author{
  Christian Koller
}
\name{torrentR-package}
\alias{torrentR-package}
\docType{package}
\title{
  A package for handling data from Ion Torrent Systems
}
\description{
  This package is for exploratory analysis and visualization of data from Ion Torrent Systems.  Functionality is provided
  for loading processed signal data (from "1.wells" files) and for perfoming various subsequent anlayses and visualizations
  based.  There are plans to expand the package in the future to also handle data from the raw image files from which the
  signals are derived.
}
\details{
  \tabular{ll}{
    Package: \tab torrentR\cr
    Type: \tab Package\cr
    Version: \tab 0.4.1\cr
    Date: \tab 2011-06-20\cr
    License: \tab GPL-2\cr
    LazyLoad: \tab yes\cr
  }
}
\author{
  Simon Cawley

  Maintainer: <scawley@iontorrent.com>
}
\keyword{ package }
\name{SimulateCAFIE}
\alias{SimulateCAFIE}
\title{
  Simulate incorporation signal as a funciton of sequence and phasing parameters
}
\description{
  Given a DNA seqeunce and parameters describing the phasing process, returns
  simulated incorporate signal and associated quantities.  The function calls one
  of two different underlying models and some of the inputs are relevant for one
  and not the other - as indicated in the details below.
}
\usage{
  SimulateCAFIE(
    seq,
    flowOrder,
    cf,
    ie,
    dr,
    nflows,
    hpScale      = 1,
    simModel     = c("treePhaserSim","CafieSolver","PhaseSim"),
    hpSignal     = 0:7,
    sigMult      = 1,
    conc         = diag(4),
    maxAdvances  = 2,
    droopType    = c("ONLY_WHEN_INCORPORATING","EVERY_FLOW"),
    extraTaps    = 0,
    getStates    = 0,
    diagonalStates = 0,
    RecalModelFile="",
    RecalModelThreshold=4,
    xval=NA,
    yval=NA
  )
}
\arguments{
  \item{seq}{
    The DNA sequence for which to simulate.
  }
  \item{flowOrder}{
    The flow cycle - for example "TACG".
  }
  \item{cf,ie,dr}{
    Estimates for cf, ie and dr.  For simModel="CafieSolver" these must be be scalars.  For
    simModel="PhaseSim" they can be scalars or vectors.  If vectors, then values will be
    cycled over flows.
  }
  \item{nflows}{
    The number of flows to model.
  }
  \item{hpScale}{
    HpScaling factor - incorporation signals for an HP of length h will be modeled as h*hpScale^(h-1).  Can be of length 1 or 4, in the
    case of the latter it is interpreted as a vector of per-nuc values in the order A,C,G,T.
  }
  \item{simModel}{
    The simulation engine to use - can be either "CafieSolver" or "PhaseSim" or "treePhaserSim".  In the case of the
    first, simulation is performed by CafieSolver::SimulateCAFIE(), for the second PhaseSim::simulate() is used.
  }
  \item{hpSignal}{
    The homopolymer response - only used when simModel="CafieSolver".  Can be used to model
    nonlinear homopolymer response.  The first element gives the expected response for a 0mer,
    the second for a 1mer, and so on.  The number of elements should be equal to MAX_MER as
    defined in CafieSolver.h (or a warning will be delivered).
  }
  \item{sigMult}{
    Value by which to multiply all incorporation signals - only used when simModel="CafieSolver".
  }
  \item{conc}{
    The 4x4 nucleotide concentration matrix to use.  Column and row order is ACGT.  The value in
    row i and colum j is the amount of nucleotide j that is present when flowing nucleotide i.
    The default is to use the identity matrix.
  }
  \item{maxAdvances}{
    The maximum number of homopolymer stretches that can be extended in a single flow.
  }
  \item{droopType}{
    The droop model used - can be either "ONLY_WHEN_INCORPORATING" (the default) or "EVERY_FLOW".
  }
  \item{extraTaps}{
    Controls the amount of extra flows to apply after each nuc flow.  The idea is to model situations where
    extra flows are applied to try drive to complete extension, though signal isn't actually collected on these
    flows.
  }
    \item{getStates}{
    Advanced simulation that also logs the state of each incorporating homopolymer.
  }
    \item{diagonalStates}{
    Switch to enable a diagonal state progression model (only with treephaser simulator).
  }
  \item{RecalModelFile}{
    Filename of a HP recalibration file to be loaded. (only with treephaser simulator)
  }
  \item{RecalModelThreshold}{
    Lower (inclusive) threshold for model HP recalibration. (default 4)
  }
  \item{xval}{
    x coordinates of the wells. Required when RecalModelFile is provided.
  }
  \item{yval}{
    y coordinates of the wells. Required when RecalModelFile is provided.
  }
}
\value{
  The return value is a list, some of whose element depend on which simulation model is being used.
  \item{sig}{
    The simulated incorporation signal
  }
  \item{hpWeight}{
    Only returned when simModel="PhaseSim".  Returns a matrix recording the proportion of live templates in
    each state for each flow of the simulation.  One row per flow and one column for each state, where the
    first state is the null state (no incorporation) and the last state is full template incorporation.
  }
  \item{droopWeight}{
    Only returned when simModel="PhaseSim".  Identifies the proporation of drooped templates at each flow.
  }
}
\seealso{
  \code{\link{phaseFit}}, \code{\link{phaseSolve}},
}
\author{
  Stuart Davidson and Simon Cawley
}
System Administration and Configuration
=======================================

In order to fully enable all Torrent PC Analysis Suite features, some
system configuration may be necessary.

Ports
-----

The following ports must allow incoming and outgoing traffic from the 
database server.

* **22**: Secure Shell (SSH)
* **80**: HTTP
* **443**: HTTPS/SSL
* **6444**: Sun Grid Engine Master
* **6445**: Sun Grid Engine Execution Daemon

Permissions
-----------

* All files and folders in ``/opt/ion/iondb/`` should be owned by ``www-data``.
* Folder in the analysis output directory currently need to have permissions
  set to 755. This is so that Sun Grid Engine jobs can write to the
  directories (as ``ion``) while simultaneously allowing Apache to serve
  the contents of those directories (as ``www-data``).

Security
--------

In order to restrict access to the database, you will need to set up a
combination of HTTP Basic Authentication and HTTPS.


Passwords
---------

By convention, the UI uses the username/password combination
*ionadmin* / *ionadmin*. These credentials will get you through basic
authentication at collaborator sites (the password is different internally),
and will also give you access to the `administrative interface
<https://analysis.iontorrents.com/admin/>`_ at all sites.

The Linux usernames on collaborator's machines are not quite as standardized.
You can attempt to SSH in with username *ion* and the Standard Password,
or try *ionguest* / *ionguest*... automodule:: rundb.ajax
   :members:
   :undoc-members:.. automodule:: crawler
   :members:
   :undoc-members:User Interface Internals
========================

The Torrent PC Analysis Suite user interface mixes server-side
`Django <http://djangoproject.com/>`_
templating with client side JavaScript, powered by `jQuery
<http://jquery.com/>`_. The entire interface is styled with CSS found in
``media/stylesheet.css``. Furthermore, some styling is applied directly with
JavaScript. As one might guess, these components do not always interact in
the most transparent ways.

The UI's tooltips provide a good example. Here is the entire process required
to pop up a tooltip window when the user hovers over a tooltip-enabled piece
of text. Let's consider the "search by experiment" tooltip found on the
`experiments <https://analysis.iontorrents.com/rundb/>`_ page.

#. When the web server starts, it imports the ``iondb.rundb.tooltips`` module,
   which in turn loads the file ``iondb/rundb/tooltips.txt``, and parses
   the tooltip descriptions it contains. One of these descriptions is titled
   "Search By Experiment" and has the key "exp_searchbyexp."
#. When a user requests the experiments page, Django renders a template,
   specifically ``iondb/templates/rundb/ion_experiment_new.html``. This
   template contains a custom `Django tag
   <http://docs.djangoproject.com/en/dev/howto/custom-template-tags/>`_
   called ``tooltip``, which is loaded with ``{% load embeddedhelp %}``.
   The tag ``{% tooltip exp_searchbyexp %}`` generates a snippet of
   html::
   
     <span class="tooltip">/rundb/tooltip/exp_searchbyexp</span>

#. The CSS class "tooltip" makes this ``span`` invisible (``display:none``).
#. When the browser has finished loading the page, it fires a series of
   `jQuery callbacks <http://api.jquery.com/jQuery/#jQuery3>`_ that
   add additional styling to the page and attach event handlers to certain
   nodes in the `DOM <http://en.wikipedia.org/wiki/Document_Object_Model>`_.
   One of these callbacks handles tooltips. It finds all `DOM nodes with
   class "tooltip" <http://api.jquery.com/jQuery/#jQuery1>`_
   (``$(".tooltip")``).
   It then looks at their direct parents in the DOM tree
   (``$(".tooltip").parent())``), and applies the CSS class "tooltip_parent" to
   these nodes. As a result, the text "Search by Experiment" is underlined with
   a fine dotted black line.
#. Another jQuery callback attaches a ``hover`` event to the tooltip parent
   nodes. The hover event, when triggered, calls the function
   ``retrieveTooltip``.
#. When the user hovers (mouses over) the underlined "Search by Experiment"
   text, the ``retrieveTooltip`` function is called. It first applies the
   "tooltip_highlighted" CSS class to the "Search by Experiment" text. Then
   it fires off an AJAX request to "/rundb/tooltip/exp_searchbyexp", the
   text contained in the "exp_searchbyexp" tooltip span.
#. The web server receives this request, which is processed by the view
   ``iondb.rundb.ajax.tooltip``. The last part of the URL, "exp_searchbyexp"
   tells the ``tooltip`` view the key of the tooltip to return. The ``tooltip``
   view then looks for "exp_searchbyexp" by calling
   ``iondb.rundb.tooltips.tip``. The function returns a Python object which
   can be serialized to `JSON <http://json.org>`_. Usually, this will
   be a dictionary containing a "title" field and a "text" field. The
   ``tooltip`` view serializes this dictionary to a UTF-8 encoded JSON string
   and returns it as a response to the original AJAX request.
#. Upon receiving the response, the browser calls a JavaScript closure
   generated by the function ``tooltipCbFactory``. This closure creates
   a `jQuery UI dialog box <http://jqueryui.com/demos/dialog/>`_. The
   closure then sets the dialog's title and content according to the
   "title" and "text" fields it received in the AJAX response. Finally,
   it displays the dialog to the user slightly below the "Search By Experiment"
   text.

jQuery Callbacks
----------------

These JavaScript functions are called by ``$(document).ready()``. They appear
in ``iondb/media/scripts.js``. These functions follow the naming convention
``prep_*``.

.. function:: prep_graphable()

   Generates an ionogram from a test fragment template. For each element
   with class ``graphable``, the function looks for this pattern in the DOM: ::

     <tr class="row1 graphable" id="1_row">...</tr>
     <tr id="1">
       <td class="sequence">TGTGA...</td>
       ...
       <div  id="1_graph" class="graph_holder"></div>
     </tr>
   
   Notice how the first ``tr`` has ID "1_row" and the next table row has
   ID "1". The first row ID must start with "<row number>_" in order to 
   be matched with the sequence in the net row.
   
   ``prep_graphable`` calls :func:`graph_template` to generate an ionogram,
   with the ``sequence`` argument extracted from the node with class
   "sequence".
   :func:`graph_template` then inserts the ionogram image into the node with
   ID "<row number>_graph".

.. function:: prep_tooltip()

   First, adds the class "tooltip_parent" to the parent nodes of everything
   matched
   by ``$(".tooltip")``. This adds the dotted underline to all tooltip text.
   
   Next, ``prep_tooltip`` binds the ``retrieveTooltip`` function to the
   "hover" event. When the event is triggered, ``retrieveTooltip`` begins
   the process of getting the tooltip from the web server.

.. function:: prep_tooltip_summary()

   Binds the function ``gatherTooltips`` on click to DOM nodes with class
   "tooltip_summary". When the event is triggered, ``gatherTooltips`` retrieves
   all tooltips on the page from the web server and displays them in a single
   large `modal <http://jqueryui.com/demos/dialog/#option-modal>`_ dialog box.

.. function:: prep_controlform()

   Binds the function ``submitControlForm`` on change to all ``input`` and
   ``select`` elements within a DOM subtree of the form::

     <form id="#control_form" ... >
       <table><tbody><tr><td>
         ...
       </td></tr></tbody></table>
     </form>

   Whenever an ``input`` or ``select`` is modified, the control form is
   submitted, and the list of database objects (either ``models.Experiment``
   or ``models.Results`` objects) is updated.

.. function:: prep_tabs()

   Binds a hover event to all DOM nodes with class "tabtext". When triggered,
   the event changes the node's background color from gray ("#cccccc") to
   white ("#ffffff"). The currently selected tab does not change color.

.. function:: prep_tab_corners()

   Applies two jQuery UI styles to DOM nodes with class "tabtext":
   `ui-corner-tl <http://jqueryui.com/docs/Theming/API#Corner_Radius_helpers>`_
   and `ui-corner-tr
   <http://jqueryui.com/docs/Theming/API#Corner_Radius_helpers>`_.

.. function:: prep_star()
   
   Binds the function ``star`` on click to ``input`` elements within DOM nodes
   of class "star_td". The ``start`` function sends an AJAX request to the
   web server, instructing Django to set or unset the ``star`` field of
   a ``models.Experiment`` object.

.. function:: prep_icon_toggling()
   
   Binds a click event to all DOM nodes of the form::

     <span class="icon_link">
       <span class="__icon_1">class_name_1</span>
       <span class="__icon_2">class_name_2</span>
     </span>

   Although the node types must not necessarily be ``span``. When the click
   event is triggered, the node with class "icon_link" removes
   the class contained in the text of the "__icon_1" node and adds that from the
   text of the "__icon_2" node.

   For example::
   
     <!-- before click -->
     <span class="icon_link class_name_1">
       <span class="__icon_1">class_name_1</span>
       <span class="__icon_2">class_name_2</span>
     </span>

     <!-- after click -->
     <span class="icon_link class_name_2">
       <span class="__icon_1">class_name_1</span>
       <span class="__icon_2">class_name_2</span>
     </span>

.. function:: prep_icon_effects()

   Adds highlighting on mouseover to DOM nodes with class "icon_link".

.. function:: prep_sorting_text()

   Bind a click event to DOM nodes matching
   ``$(".sortables > th > .sortheading")``. When triggered the event calls
   ``setSorting`` with the node's parent element as the argument. This allows
   for toggling ascending/descending sorting based on the clicked node (see
   `sorting` for more details).

.. function:: prep_sorting_buttons()

   Add up-arrow and down-arrow buttons for sorting to sortable elements
   (roughly, those nodes that are matched by
   ``$(".sortables > th > .sortheading")``.

.. function:: prep_centering_ie6()
   
   Add the class "centered_ie6" to nodes with class "centered" if the 
   browser does not support "box model" rendering.

.. function:: prep_centering_width()
   
   Bind an on resize event to the browser window. When the window is resized,
   the ``div`` node with class "all", which wraps the entire UI, is resized
   as well. This seems to be a nice cross-browser solution to maintaining a
   fixed-width margin on both the left and right sides of the screen.

JavaScript Functions
--------------------

.. data:: MOUSEOVER_DELAY

   Number of milliseconds between the time the user hovers over a tooltip
   and the appearance of a tooltip dialog box. Currently set to 750.

Visibility Controls
^^^^^^^^^^^^^^^^^^^

Functions for showing/hiding DOM nodes.

.. function:: toggleTr(id)

   Perform the jQuery `toggle effect <http://api.jquery.com/toggle/>`_
   on the element with ID ``id + "_holder"``::

      <a href="javascript:toggleTr('123')">
      	 Click me to reveal the hidden DIV below.
      </a>
      <div id="123_holder">Click the link above to reveal me.</div>

.. function:: toggleAdvanced(id)

   Toggle (see :func:`toggleTr`) the node with ID ``id``.

.. function:: clickExpand(id)

   Programmatically click all anchors with class "icon_link" beneath the
   node with ID ``id``

Starring
^^^^^^^^

.. function:: star()

   Dispatch an AJAX request to "star" a ``models.Experiment`` database
   record::

      <input type="checkbox" id="star_123" onchange="star()"/>
    
   Clicking this checkbox will "star" the star the ``models.Experiment`` record
   with primary key 123.

Job Control and Termination
^^^^^^^^^^^^^^^^^^^^^^^^^^^

Functions for terminating jobs. These methods are used in the "Jobs" page,
which is rendered by :func:`rundb.views.current_jobs`.

.. function:: do_control(url)

   Make an AJAX request using jQuery's `getJSON 
   <http://api.jquery.com/jQuery.getJSON/>`_ to ``url``. Registers a callback
   with `getJSON` that displays the result of the AJAX call (either 
   "termination succeded" or "termination failed") to the user.

.. function:: build_control_dialogue(url,name)

   Display a jQuery `dialog box <http://jqueryui.com/demos/dialog>` which
   presents the user with the option to terminate a job. If the user
   clicks "Terminate", the function calls :func:`do_control` with argument
   ``url``.

TF Ionogram Creation
^^^^^^^^^^^^^^^^^^^^

.. function:: graph_template(sequence,graph_node)

   Convert nucleotide sequence ``sequence`` into flowspaces using floworder
   "TACG", then create an ionogram using `jQuery Google Charts API
   <http://www.maxb.net/scripts/jgcharts/include/demo/>`_.
   :func:`graph_template` then inserts the ionogram image into the
   DOM node ``graph_node``.

Searching and Sorting Functions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. function:: submitControlForm()

   Submit the control form, generating a new set of results. This function
   is called by several sorting-related functions.

.. function:: setSorting(sortkey)

   Set the value of the :ref:`sortfield <sortfield>` to ``sortkey``. If the
   sortfield is already set to ``sortkey``, reverse the sort order.
   For example, if ``sortkey == "experiment"`` and the sortfield contains
   ``"project"`` then :func:`setSorting` will set the the sortfield to
   ``"-project"``.

   After setting the sortfield, :func:`setSorting` calls
   :func:`submitControlForm`.

.. function:: getSortKey(sortable_th)

   Extract the :ref:`sortkey <sortkey>` from the jQuery-wrapped ``th`` node
   ``sortable_th``.

.. function:: selectedSortField()

   Return the :ref:`sortkey <sortkey>` corresponding to the current field
   for sorting,
   extracted from the :ref:`sortfield <sortfield>`. For example, if
   the sortfield
   contains either ``"project"`` or ``"-project"``, :func:`selectedSortField`
   returns ``"project"``.

.. function:: sorterIsSelected(sortkey)

   Return ``true`` if ``sortkey`` is the current field selected for sorting.

.. function:: sortClickCbFactory(rev,sortkey,sfelement)

   Return a closure which sets the value of ``sfelement`` to ``'-' + sortkey``
   or ``sortkey`` depending on whether ``rev`` is ``true`` or ``false``,
   respectively.

   Normally, ``sfelement`` should be the :ref:`sortfield <sortfield>`.

.. function:: addSortIcon()

   Creates two `icons <http://jqueryui.com/themeroller/>`_, an up arrow
   (``ui-icon-arrowthick-1-n``) and a down arrow (``ui-icon-arrowthick-1-s``).
   :func:`addSortIcon` then appends the icons to the node with class 
   "sortheading" under node ``this``.

.. function:: docHasSortable()
   
   Returns ``true`` if the page contains a :ref:`sortfield <sortfield>`.

Tooltip Functions
^^^^^^^^^^^^^^^^^

.. data:: HOVER_TIMEOUTS

   Dictionary mapping timeout ID's returned by `setTimeout
   <http://www.w3schools.com/js/js_timing.asp>` to URL's which identify the
   tooltip to be displayed once the timeout has expired. Keeping around
   these mappings makes it easy to cancel a timeout in the event that the
   user mouses off of a tooltip before the timeout has passed.

.. data:: TOOLTIP_KEYS

   A cache of previously loaded tooltip data.

.. function:: tooltipCbFactory(ele,url)

   Return a closure takes a JSON object returned by jQuery `getJSON()
   <http://api.jquery.com/jQuery.getJSON/>`_ and displays the tooltip
   information in that object using a dialog box. The dialog box appears
   slightly below DOM node ``ele``.

.. function:: tooltipClose()

   Remove the "tooltip_highlighted" class from ``this``, cancel any
   timeouts associated with the tooltip under ``this``, and hide
   the tooltip dialog box if it is showing.

.. function:: extractUrl(ele)

   Return the tooltip URL contained in the text of ``ele``. This URL can
   be used to:
   
   * Uniquely identify a tooltip
   * Retrieve data from the web server (via jQuery `getJSON 
     <http://api.jquery.com/jQuery.getJSON/>`_) to display the tooltip.

   :func:`extractUrl` extracts the tooltip URL from the following structure::

     <span id="the_ele_passed_in>
         ...
	 <span class=".tooltip">/the/url/to/extract</span>
     </span>	   
    
.. function:: _tt(url,cb)

   Checks :data:`TOOLTIP_KEYS` to see if the tooltip for ``url`` has already
   been retrieved. If so, call ``cb`` on the the tooltip data. Otherwise,
   use jQuery `getJSON() <http://api.jquery.com/jQuery.getJSON/>`_ to
   retrieve the tooltip data, using ``cb`` as the callback argument to
   `getJSON()`.

.. function:: retrieveTooltip()

   First, adds the class "tooltip_highlighted" to ``this``. Next, calls 
   `setTimeout <http://www.w3schools.com/js/js_timing.asp>`_, which begins
   the process of retrieving and displaying the tooltip for ``this``.

.. function:: displaySummary(d,keys)

   Given the data for all tooltips on a page (``d``), and the URL's for each
   tooltip in ``keys``, pop up a dialog box that displays the title and
   text of every tooltip on the page.

.. function:: retrieveAllTooltips(tips)

   For each element in ``tips`` retrieve the associated tooltip data and
   store the data in a dictionary keyed by the tooltip's URL. When all
   tooltips have been entered into the dictionary, call :func:`displaySummary`
   with the dictionary and list of URL's as arguments.

.. gatherTooltips()
   
   Find all tooltips on the page with ``$(".tooltip).parent()`` and
   call :func:`retrieveAllTooltips` on the result.
   

Sorting
-------

The sortable columns in the "Experiments" and "Reports" pages rely on elements
from all levels of the UI stack.

Making A Column Sortable
^^^^^^^^^^^^^^^^^^^^^^^^

The UI's sorting mechanism assumes that you are attempting to sort data laid
out in an HTML table. There should be a ``th`` elements at the top of each
column. The top row (``tr``) of the column must have the class "sortables".
Each ``th`` at the top of a sortable column must contain at least the following
two nodes:

* A node with class "sortheading"
* A node with class "sortkey"

The structure of the table should look something like this::

  <table><thead>
    <tr class="sortables">
      ...
      <th>
        <div class="sortheading">Name of Column</div>
	<div class="sortkey">corresponding_model_field</div>
	<!-- additional nodes allowed here -->
      </th>
     </tr>
   </thead>...</table>

.. _sortkey:

The "sortheading" node should contain the name of the column - the actual text
that will appear on the rendered page. The "sortkey" node will be hidden. It
should contain the name of the database field corresponding that will actually
be used for sorting. For example, when sorting ``model.Experiment`` objects::

   <div class="sortheading">PGM</div>
   <div class="sortkey">pgmname</div>

could be used to sort by the ``pgmName`` setting. Note that the "sortkey" text
is case-insensitive.

The :func:`prep_sorting_text` and :func:`prep_sorting_buttons` jQuery callbacks
will take care of adding the necessary buttons and event handlers to make the
column sortable.

Enabling Sorting On A Page
^^^^^^^^^^^^^^^^^^^^^^^^^^

The column to sort on, and the sort order (ascending vs. descending) is passed
to the web server with a hidden form field. The page containing the sortable
table should include a ``form`` with ID "control_form". That is, it should
contain something like::

   <form id="control_form" ...>
      ...
   </form>

.. _sortfield:

The view that renders the page's template (for example
:func:`rundb.views.experiment`) should create an instance of
``forms.SortForm``. The template should then render the ``SortForm``'s
``sortfield`` within the control form::

   <form id="control_form" ...>
     ...
     <!-- other forms and fields -->
     ...
     {{sortform.sortfield}}
   </form>

This example assumes that the ``SortForm`` instance is included in the
template context as "sortform". The field inserted into the form is styled
so as to be invisible -- it uses input type "hidden".

In order to actually sort information, the view will need to process the data
found in the sortfield. The :mod:`rundb.views` module provides an easy way to
sort (and search) data for a particular model. Use the
:func:`rundb.views.search_and_sort`
function, with a ``SortForm`` instance as one of the arguments, to return
a sorted queryset.

Sorting Internals
^^^^^^^^^^^^^^^^^

The :func:`prep_sorting_text` and :func:`prep_sorting_buttons` significantly
modify the page's DOM tree in order to add UI elements necessary for sorting.

:func:`prep_sorting_text` binds a click event to sortable column titles. When
the event is triggered, it calls :func:`setSorting` with the column's *sortkey*
as the argument. The *sortkey* is the content of the node with class "sortkey"
in the ``th`` element holding the column name::

   ...
   <th>
     <div class="sortheading">The Column Name</div>
     <div class="sortkey">the_sortkey</div>
   </th>
   ...

:func:`setSorting` sets the value of the the hidden input with class
"sortfield" living in the control form. It then calls
:func:`submitControlForm`, which submits the form.

:func:`prep_sorting_buttons` calls :func:`addSortIcon` on all nodes matching
``$(".sortables > th")``.

Searching
---------

As with sorting, searching assumes the presence of a form with ID
"control_form". A search-enabled view (such as :func:`rundb.views.experiment`)
should create an instance of :class:`rundb.forms.SearchForm`.

The template should render the ``SearchForm``'s ``searchterms`` field::

   <form id="control_form" ...>
     <table>
       ...
         <tr>
	 ...
	   <td>{{searchform.searchterms}}</td>
	 ...
	 </tr>
     </table> 
   </form>

This assumes that the ``SearchForm`` instance is included in the template
context as ``searchform``. The ``searchterms`` field has class "searchbox",
which appears in ``iondb/media/stylesheet.css``.

In the view, the ``SearchForm`` instance should be passed to
:func:`rundb.views.search_and_sort` for the purposes of searching the
queryset.

CSS
---

The UI's styling is probably its least organized component. There is no real
naming convention or grouping to all the styles in
``iondb/media/stylesheet.css``. That being said, it should be a fairly
straightforward job to refactor all the styling (see the :ref:`future-work`
section).

In addition, styles used in the Django templates mix those defined in
``stylesheet.css`` and those defined by the UI's jQuery theme. Anything
that matches the wildcard pattern ``ui-*`` is generally a jQuery style.

Quick Tips and Gotchas
----------------------

* If you add a tooltip to a template and get an "invalid tag" Django exception,
  you need to add the tag ``{% load embeddedhelp %}`` somewhere in that
  template.
* Much of the CSS for buttons is applied by JavaScript. If you're trying to
  change the look of an interactive component and modifying
  ``iondb/media/stylesheet.css``
  isn't having an effect, you might want to check ``iondb/media/scripts.js``.
  While handling some styling in JavaScript made it easier to play with look
  and feel, this should probably be factored out into a stylesheet.

.. _future-work:

Future Work
-----------

* Refactor ``iondb/media/scripts.js`` into several files, with the files
  organized by functionality. For example, we might want files such as
  ``iondb/media/sorting.js`` and ``iondb/media/searching.js``.
* Add namespacing to ``iondb/media/scripts.js``.
* Organize ``iondb/media/stylesheet.css`` and perhaps rename styles in order
  to achieve a consistent naming convention. Break the file up into subfiles.
* `Roll a custom jQuery theme <http://jqueryui.com/themeroller/>`_.
* Replace the two-button sorting controls with a single, toggling button... automodule:: rundb.models
   :members:
   :undoc-members:.. automodule:: serve
   :members:
   :undoc-members:.. Torrent PC Analysis Suite documentation master file, created by sphinx-quickstart on Tue Jan 19 15:50:59 2010.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Torrent PC Analysis Suite Documentation
=======================================

.. toctree::
   :maxdepth: 2

   serve
   crawler
   views
   ajax
   models
   ui
   sysadmin

Emergency Contact Info
======================

You can reach Jeremy by e-mail at ``jhoon at iontorrents`` or on his
mobile at (203) 215-5031.

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

.. automodule:: rundb.views
   :members:
   :undoc-members: