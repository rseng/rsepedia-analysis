DateTime form input field for ExtJS.

Build
-----

Build requires:
* SenchaCmd, download from [Sencha](http://www.sencha.com/products/sencha-cmd/download/) and install.
* karma, for running tests, install with `sudo apt-get install nodejs lcov` and `sudo npm install -g karma-cli`.
* karma plugins, `cd packages/datetime/;npm install`.
* jsduck, for documentation, install with `sudo gem install jsduck`.

Build with `sencha package build` inside `packages/datetime/` folder.

To make pkg available run `sencha package add ../../build/datetime/datetime.pkg` and make `<SenchaCmd Installation Prefix>/repo/pkgs` folder available online.

Installation and usage
----------------------

To use with a SenchaCmd based project:

1. `sencha repo add -address <to be announced> -name NLeSC
2. Add 'datetime' to `requires` array in `app.json`.
3. Use component in your application.
4. `sencha app refresh` do download package.

To use with a non-SenchaCmd project:

1. unzip `datetime.pkg`
2. Add `build/datetime.js` into html page.

Documentation
-------------

Documentation can be generated with `ant docs`.
To make inline examples work, make a copy of the ExtJS SDK directory as `docs/extjs-build`.
Place the docs folder online and open `docs/index.html` in a webbrowser.

Copyrights & Disclaimers
------------------------

extjs-datetime is copyrighted by the Netherlands eScience Center and releases under
the Apache License, Version 2.0.

See <http://www.esciencecenter.nl> for more information on the Netherlands
eScience Center.

See the "LICENSE" file for more information.
# datetime - Read Me

# datetime/src

This folder contains source code that will automatically be added to the classpath when
the package is used.
# datetime/sass

This folder contains SASS files of various kinds, organized in sub-folders:

    datetime/sass/etc
    datetime/sass/src
    datetime/sass/var
# datetime/sass/src

This folder contains SASS sources that mimic the component-class hierarchy. These files
are gathered in to a build of the CSS based on classes that are used by the build.
# datetime/sass/var

This folder contains variable declaration files named by their component class.
# datetime/sass/etc

This folder contains miscellaneous SASS files. Unlike `"datetime/sass/etc"`, these files
need to be used explicitly.
# datetime/overrides

This folder contains overrides which will automatically be required by package users.
# datetime/resources

This folder contains static resources (typically an `"images"` folder as well).
