
# About

PotreeConverter generates an octree LOD structure for streaming and real-time rendering of massive point clouds. The results can be viewed in web browsers with [Potree](https://github.com/potree/potree) or as a desktop application with [PotreeDesktop](https://github.com/potree/PotreeDesktop). 

Version 2.0 is a complete rewrite with following differences over the previous version 1.7:

* About 10 to 50 times faster than PotreeConverter 1.7 on SSDs.
* Produces a total of 3 files instead of thousands to tens of millions of files. The reduction of the number of files improves file system operations such as copy, delete and upload to servers from hours and days to seconds and minutes. 
* Better support for standard LAS attributes and arbitrary extra attributes. Full support (e.g. int64 and uint64) in development.
* Optional compression is not yet available in the new converter but on the roadmap for a future update.

Altough the converter made a major step to version 2.0, the format it produces is also supported by Potree 1.7. The Potree viewer is scheduled to make the major step to version 2.0 in 2021, with a rewrite in WebGPU. 

# Publications

* [Potree: Rendering Large Point Clouds in Web Browsers](https://www.cg.tuwien.ac.at/research/publications/2016/SCHUETZ-2016-POT/SCHUETZ-2016-POT-thesis.pdf)
* [Fast Out-of-Core Octree Generation for Massive Point Clouds](https://www.cg.tuwien.ac.at/research/publications/2020/SCHUETZ-2020-MPC/), _Schütz M., Ohrhallinger S., Wimmer M._

# Getting Started

1. Download windows binaries or
    * Download source code
	* Install [CMake](https://cmake.org/) 3.16 or later
	* Create and jump into folder "build"
	    ```
	    mkdir build
	    cd build
	    ```
	* run 
	    ```
	    cmake ../
	    ```
	* On linux, run: ```make```
	* On windows, open Visual Studio 2019 Project ./Converter/Converter.sln and compile it in release mode
2. run ```PotreeConverter.exe <input> -o <outputDir>```
    * Optionally specify the sampling strategy:
	* Poisson-disk sampling (default): ```PotreeConverter.exe <input> -o <outputDir> -m poisson```
	* Random sampling: ```PotreeConverter.exe <input> -o <outputDir> -m random```

In Potree, modify one of the examples with following load command:

```javascript
let url = "../pointclouds/D/temp/test/metadata.json";
Potree.loadPointCloud(url).then(e => {
	let pointcloud = e.pointcloud;
	let material = pointcloud.material;

	material.activeAttributeName = "rgba";
	material.minSize = 2;
	material.pointSizeType = Potree.PointSizeType.ADAPTIVE;

	viewer.scene.addPointCloud(pointcloud);
	viewer.fitToScreen();
});

```

# Alternatives

PotreeConverter 2.0 produces a very different format than previous iterations. If you find issues, you can still try previous converters or alternatives:

<table>
	<tr>
		<th></th>
		<th>PotreeConverter 2.0</th>
		<th><a href="https://github.com/potree/PotreeConverter/releases/tag/1.7">PotreeConverter 1.7</a></th>
		<th><a href="https://entwine.io/">Entwine</a></th>
	</tr>
	<tr>
		<th>license</th>
		<td>
			free, BSD 2-clause
		</td>
		<td>
			free, BSD 2-clause
		</td>
		<td>
			free, LGPL
		</td>
	</tr>
	<tr>
		<th>#generated files</th>
		<td>
			3 files total
		</td>
		<td>
			1 per node
		</td>
		<td>
			1 per node
		</td>
	</tr>
	<tr>
		<th>compression</th>
		<td>
			none (TODO)
		</td>
		<td>
			LAZ (optional)
		</td>
		<td>
			LAZ
		</td>
	</tr>
</table>

Performance comparison (Ryzen 2700, NVMe SSD):

![](./docs/images/performance_chart.png)

# License 

PotreeConverter is available under the [BSD 2-clause license](./LICENSE).<p align="center"><img src="https://brotli.org/brotli.svg" alt="Brotli" width="64"></p>

# SECURITY NOTE

Please consider updating brotli to version 1.0.9 (latest).

Version 1.0.9 contains a fix to "integer overflow" problem. This happens when "one-shot" decoding API is used (or input chunk for streaming API is not limited), input size (chunk size) is larger than 2GiB, and input contains uncompressed blocks. After the overflow happens, `memcpy` is invoked with a gigantic `num` value, that will likely cause the crash.

### Introduction

Brotli is a generic-purpose lossless compression algorithm that compresses data
using a combination of a modern variant of the LZ77 algorithm, Huffman coding
and 2nd order context modeling, with a compression ratio comparable to the best
currently available general-purpose compression methods. It is similar in speed
with deflate but offers more dense compression.

The specification of the Brotli Compressed Data Format is defined in [RFC 7932](https://tools.ietf.org/html/rfc7932).

Brotli is open-sourced under the MIT License, see the LICENSE file.

Brotli mailing list:
https://groups.google.com/forum/#!forum/brotli

[![TravisCI Build Status](https://travis-ci.org/google/brotli.svg?branch=master)](https://travis-ci.org/google/brotli)
[![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/github/google/brotli?branch=master&svg=true)](https://ci.appveyor.com/project/szabadka/brotli)
[![Fuzzing Status](https://oss-fuzz-build-logs.storage.googleapis.com/badges/brotli.svg)](https://oss-fuzz-build-logs.storage.googleapis.com/index.html#brotli)

### Build instructions

#### Vcpkg

You can download and install brotli using the [vcpkg](https://github.com/Microsoft/vcpkg/) dependency manager:

    git clone https://github.com/Microsoft/vcpkg.git
    cd vcpkg
    ./bootstrap-vcpkg.sh
    ./vcpkg integrate install
    vcpkg install brotli

The brotli port in vcpkg is kept up to date by Microsoft team members and community contributors. If the version is out of date, please [create an issue or pull request](https://github.com/Microsoft/vcpkg) on the vcpkg repository.

#### Autotools-style CMake

[configure-cmake](https://github.com/nemequ/configure-cmake) is an
autotools-style configure script for CMake-based projects (not supported on Windows).

The basic commands to build, test and install brotli are:

    $ mkdir out && cd out
    $ ../configure-cmake
    $ make
    $ make test
    $ make install

By default, debug binaries are built. To generate "release" `Makefile` specify `--disable-debug` option to `configure-cmake`.

#### Bazel

See [Bazel](http://www.bazel.build/)

#### CMake

The basic commands to build and install brotli are:

    $ mkdir out && cd out
    $ cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=./installed ..
    $ cmake --build . --config Release --target install

You can use other [CMake](https://cmake.org/) configuration.

#### Premake5

See [Premake5](https://premake.github.io/)

#### Python

To install the latest release of the Python module, run the following:

    $ pip install brotli

To install the tip-of-the-tree version, run:

    $ pip install --upgrade git+https://github.com/google/brotli

See the [Python readme](python/README.md) for more details on installing
from source, development, and testing.

### Benchmarks
* [Squash Compression Benchmark](https://quixdb.github.io/squash-benchmark/) / [Unstable Squash Compression Benchmark](https://quixdb.github.io/squash-benchmark/unstable/)
* [Large Text Compression Benchmark](http://mattmahoney.net/dc/text.html)
* [Lzturbo Benchmark](https://sites.google.com/site/powturbo/home/benchmark)

### Related projects
> **Disclaimer:** Brotli authors take no responsibility for the third party projects mentioned in this section.

Independent [decoder](https://github.com/madler/brotli) implementation by Mark Adler, based entirely on format specification.

JavaScript port of brotli [decoder](https://github.com/devongovett/brotli.js). Could be used directly via `npm install brotli`

Hand ported [decoder / encoder](https://github.com/dominikhlbg/BrotliHaxe) in haxe by Dominik Homberger. Output source code: JavaScript, PHP, Python, Java and C#

7Zip [plugin](https://github.com/mcmilk/7-Zip-Zstd)

Dart [native bindings](https://github.com/thosakwe/brotli)
Want to contribute? Great! First, read this page (including the small print at
the end).

### Before you contribute
Before we can use your code, you must sign the
[Google Individual Contributor License Agreement]
(https://cla.developers.google.com/about/google-individual)
(CLA), which you can do online. The CLA is necessary mainly because you own the
copyright to your changes, even after your contribution becomes part of our
codebase, so we need your permission to use and distribute your code. We also
need to be sure of various other things—for instance that you'll tell us if you
know that your code infringes on other people's patents. You don't have to sign
the CLA until after you've submitted your code for review and a member has
approved it, but you must do it before we can put your code into our codebase.
Before you start working on a larger contribution, you should get in touch with
us first through the issue tracker with your idea so that we can help out and
possibly guide you. Coordinating up front makes it much easier to avoid
frustration later on.

### Code reviews
All submissions, including submissions by project members, require review. We
use Github pull requests for this purpose.

### The small print
Contributions made by corporations are covered by a different agreement than
the one above, the [Software Grant and Corporate Contributor License Agreement]
(https://cla.developers.google.com/about/google-corporate).
brotli(1) -- brotli, unbrotli - compress or decompress files
================================================================

SYNOPSIS
--------

`brotli` [*OPTION|FILE*]...

`unbrotli` is equivalent to `brotli --decompress`

DESCRIPTION
-----------
`brotli` is a generic-purpose lossless compression algorithm that compresses
data using a combination of a modern variant of the **LZ77** algorithm, Huffman
coding and 2-nd order context modeling, with a compression ratio comparable to
the best currently available general-purpose compression methods. It is similar
in speed with deflate but offers more dense compression.

`brotli` command line syntax similar to `gzip (1)` and `zstd (1)`.
Unlike `gzip (1)`, source files are preserved by default. It is possible to
remove them after processing by using the `--rm` _option_.

Arguments that look like "`--name`" or "`--name=value`" are _options_. Every
_option_ has a short form "`-x`" or "`-x value`". Multiple short form _options_
could be coalesced:

* "`--decompress --stdout --suffix=.b`" works the same as
* "`-d -s -S .b`" and
* "`-dsS .b`"

`brotli` has 3 operation modes:

* default mode is compression;
* `--decompress` option activates decompression mode;
* `--test` option switches to integrity test mode; this option is equivalent to
  "`--decompress --stdout`" except that the decompressed data is discarded
  instead of being written to standard output.

Every non-option argument is a _file_ entry. If no _files_ are given or _file_
is "`-`", `brotli` reads from standard input. All arguments after "`--`" are
_file_ entries.

Unless `--stdout` or `--output` is specified, _files_ are written to a new file
whose name is derived from the source _file_ name:

* when compressing, a suffix is appended to the source filename to
  get the target filename
* when decompressing, a suffix is removed from the source filename to
  get the target filename

Default suffix is `.br`, but it could be specified with `--suffix` option.

Conflicting or duplicate _options_ are not allowed.

OPTIONS
-------

* `-#`:
    compression level (0-9); bigger values cause denser, but slower compression
* `-c`, `--stdout`:
    write on standard output
* `-d`, `--decompress`:
    decompress mode
* `-f`, `--force`:
    force output file overwrite
* `-h`, `--help`:
    display this help and exit
* `-j`, `--rm`:
    remove source file(s); `gzip (1)`-like behaviour
* `-k`, `--keep`:
    keep source file(s); `zstd (1)`-like behaviour
* `-n`, `--no-copy-stat`:
    do not copy source file(s) attributes
* `-o FILE`, `--output=FILE`
    output file; valid only if there is a single input entry
* `-q NUM`, `--quality=NUM`:
    compression level (0-11); bigger values cause denser, but slower compression
* `-t`, `--test`:
    test file integrity mode
* `-v`, `--verbose`:
    increase output verbosity
* `-w NUM`, `--lgwin=NUM`:
    set LZ77 window size (0, 10-24) (default: 22); window size is
    `(2**NUM - 16)`; 0 lets compressor decide over the optimal value; bigger
    windows size improve density; decoder might require up to window size
    memory to operate
* `-S SUF`, `--suffix=SUF`:
    output file suffix (default: `.br`)
* `-V`, `--version`:
    display version and exit
* `-Z`, `--best`:
    use best compression level (default); same as "`-q 11`"

SEE ALSO
--------

`brotli` file format is defined in
[RFC 7932](https://www.ietf.org/rfc/rfc7932.txt).

`brotli` is open-sourced under the
[MIT License](https://opensource.org/licenses/MIT).

Mailing list: https://groups.google.com/forum/#!forum/brotli

BUGS
----
Report bugs at: https://github.com/google/brotli/issues
three.js
========

[![NPM package][npm]][npm-url]
[![Build Size][build-size]][build-size-url]
[![Build Status][build-status]][build-status-url]
[![Dependencies][dependencies]][dependencies-url]
[![Dev Dependencies][dev-dependencies]][dev-dependencies-url]
[![Language Grade][lgtm]][lgtm-url]

#### JavaScript 3D library ####

The aim of the project is to create an easy to use, lightweight, 3D library with a default WebGL renderer. The library also provides Canvas 2D, SVG and CSS3D renderers in the examples.

[Examples](http://threejs.org/examples/) &mdash;
[Documentation](http://threejs.org/docs/) &mdash;
[Wiki](https://github.com/mrdoob/three.js/wiki) &mdash;
[Migrating](https://github.com/mrdoob/three.js/wiki/Migration-Guide) &mdash;
[Questions](http://stackoverflow.com/questions/tagged/three.js) &mdash;
[Forum](https://discourse.threejs.org/) &mdash;
[Gitter](https://gitter.im/mrdoob/three.js) &mdash;
[Slack](https://join.slack.com/t/threejs/shared_invite/enQtMzYxMzczODM2OTgxLTQ1YmY4YTQxOTFjNDAzYmQ4NjU2YzRhNzliY2RiNDEyYjU2MjhhODgyYWQ5Y2MyZTU3MWNkOGVmOGRhOTQzYTk)

### Usage ###

Download the [minified library](http://threejs.org/build/three.min.js) and include it in your HTML, or install and import it as a [module](http://threejs.org/docs/#manual/introduction/Import-via-modules),
Alternatively, see [how to build the library yourself](https://github.com/mrdoob/three.js/wiki/Build-instructions).

```html
<script src="js/three.min.js"></script>
```

This code creates a scene, a camera, and a geometric cube, and it adds the cube to the scene. It then creates a `WebGL` renderer for the scene and camera, and it adds that viewport to the `document.body` element. Finally, it animates the cube within the scene for the camera.

```javascript
var camera, scene, renderer;
var geometry, material, mesh;

init();
animate();

function init() {

	camera = new THREE.PerspectiveCamera( 70, window.innerWidth / window.innerHeight, 0.01, 10 );
	camera.position.z = 1;

	scene = new THREE.Scene();

	geometry = new THREE.BoxGeometry( 0.2, 0.2, 0.2 );
	material = new THREE.MeshNormalMaterial();

	mesh = new THREE.Mesh( geometry, material );
	scene.add( mesh );

	renderer = new THREE.WebGLRenderer( { antialias: true } );
	renderer.setSize( window.innerWidth, window.innerHeight );
	document.body.appendChild( renderer.domElement );

}

function animate() {

	requestAnimationFrame( animate );

	mesh.rotation.x += 0.01;
	mesh.rotation.y += 0.02;

	renderer.render( scene, camera );

}
```

If everything went well you should see [this](https://jsfiddle.net/f2Lommf5/).

### Change log ###

[Releases](https://github.com/mrdoob/three.js/releases)


[npm]: https://img.shields.io/npm/v/three.svg
[npm-url]: https://www.npmjs.com/package/three
[build-size]: https://badgen.net/bundlephobia/minzip/three
[build-size-url]: https://bundlephobia.com/result?p=three
[build-status]: https://travis-ci.org/mrdoob/three.js.svg?branch=dev
[build-status-url]: https://travis-ci.org/mrdoob/three.js
[dependencies]: https://img.shields.io/david/mrdoob/three.js.svg
[dependencies-url]: https://david-dm.org/mrdoob/three.js
[dev-dependencies]: https://img.shields.io/david/dev/mrdoob/three.js.svg
[dev-dependencies-url]: https://david-dm.org/mrdoob/three.js#info=devDependencies
[lgtm]: https://img.shields.io/lgtm/grade/javascript/g/mrdoob/three.js.svg?label=code%20quality
[lgtm-url]: https://lgtm.com/projects/g/mrdoob/three.js/
# Spectrum
## The No Hassle Colorpicker

See the demo and docs: http://bgrins.github.io/spectrum.

I wanted a colorpicker that didn't require images, and that had an API that made sense to me as a developer who has worked with color in a number of applications.  I had tried a number of existing plugins, but decided to try and make a smaller, simpler one.

I started using canvas, then switched to CSS gradients, since it turned out to be easier to manage, and provided better cross browser support.

### Basic Usage

Head over to the [docs](http://bgrins.github.io/spectrum) for more information.  There is a visual demo of the different options hosted at: http://bgrins.github.io/spectrum.

    <script src='spectrum.js'></script>
    <link rel='stylesheet' href='spectrum.css' />

    <input id='colorpicker' />

    <script>
    $("#colorpicker").spectrum({
        color: "#f00"
    });
    </script>

### npm

Spectrum is registered as package with npm.  It can be installed with:

    npm install spectrum-colorpicker

### Bower

Spectrum is registered as a package with [Bower](http://bower.io/), so it can be pulled down using:

    bower install spectrum

### Using spectrum with a CDN

CDN provided by [cdnjs](https://cdnjs.com/libraries/spectrum)

    <script src="https://cdnjs.cloudflare.com/ajax/libs/spectrum/1.8.0/spectrum.min.js"></script>
    <link rel="stylesheet" type="text/css" href="https://cdnjs.cloudflare.com/ajax/libs/spectrum/1.8.0/spectrum.min.css">

### Continuous Integration

[![Build Status](https://secure.travis-ci.org/bgrins/spectrum.png?branch=master)](http://travis-ci.org/bgrins/spectrum)

Visit https://travis-ci.org/bgrins/spectrum to view the status of the automated tests.

### Building Spectrum Locally

If you'd like to download and use the plugin, head over to http://bgrins.github.io/spectrum/ and click the 'Download Zip' button.

If you'd like to run the development version, spectrum uses Grunt to automate the testing, linting, and building.  Head over to http://gruntjs.com/getting-started for more information.  First, clone the repository, then run:

    npm install -g grunt-cli
    npm install

    # runs jshint and the unit test suite
    grunt

    # runs jshint, the unit test suite, and builds a minified version of the file.
    grunt build

### Internationalization

If you are able to translate the text in the UI to another language, please do!  You can do so by either [filing a pull request](https://github.com/bgrins/spectrum/pulls) or [opening an issue]( https://github.com/bgrins/spectrum/issues) with the translation.  The existing languages are listed at: https://github.com/bgrins/spectrum/tree/master/i18n.

For an example, see the [Dutch translation](i18n/jquery.spectrum-nl.js).
# jstree

[jsTree](http://www.jstree.com/) is jquery plugin, that provides interactive trees. It is absolutely free, [open source](https://github.com/vakata/jstree) and distributed under the MIT license.

jsTree is easily extendable, themable and configurable, it supports HTML & JSON data sources, AJAX & async callback loading.

jsTree functions properly in either box-model (content-box or border-box), can be loaded as an AMD module, and has a built in mobile theme for responsive design, that can easily be customized. It uses jQuery's event system, so binding callbacks on various events in the tree is familiar and easy.

You also get:
 * drag & drop support
 * keyboard navigation
 * inline edit, create and delete
 * tri-state checkboxes
 * fuzzy searching
 * customizable node types

_Aside from this readme you can find a lot more info on [jstree.com](http://www.jstree.com) & [the discussion group](https://groups.google.com/forum/#!forum/jstree)_.

---

<!-- MarkdownTOC depth=0 autolink=true bracket=round -->

- [Getting Started](#getting-started)
  - [Include all neccessary files](#include-all-neccessary-files)
  - [Populating a tree using HTML](#populating-a-tree-using-html)
  - [Populating a tree using an array \(or JSON\)](#populating-a-tree-using-an-array-or-json)
    - [The required JSON format](#the-required-json-format)
  - [Populating the tree using AJAX](#populating-the-tree-using-ajax)
  - [Populating the tree using AJAX and lazy loading nodes](#populating-the-tree-using-ajax-and-lazy-loading-nodes)
  - [Populating the tree using a callback function](#populating-the-tree-using-a-callback-function)
- [Working with events](#working-with-events)
- [Interacting with the tree using the API](#interacting-with-the-tree-using-the-api)
- [More on configuration](#more-on-configuration)
- [Plugins](#plugins)
  - [checkbox](#checkbox)
  - [contextmenu](#contextmenu)
  - [dnd](#dnd)
  - [massload](#massload)
  - [search](#search)
  - [sort](#sort)
  - [state](#state)
  - [types](#types)
  - [unique](#unique)
  - [wholerow](#wholerow)
  - [More plugins](#more-plugins)
- [PHP demos moved to new repository](#php-demos-moved-to-new-repository)
- [License & Contributing](#license--contributing)

<!-- /MarkdownTOC -->


---

## Getting Started

### Include all neccessary files
To get started you need 3 things in your page:
 1. jQuery (anything above 1.9.1 will work)
 2. A jstree theme (there is only one theme supplied by default)
 3. The jstree source file

```html
<script src="//cdnjs.cloudflare.com/ajax/libs/jquery/3.1.1/jquery.min.js"></script>

<link rel="stylesheet" href="//cdnjs.cloudflare.com/ajax/libs/jstree/3.3.8/themes/default/style.min.css" />
<script src="//cdnjs.cloudflare.com/ajax/libs/jstree/3.3.8/jstree.min.js"></script>
```

_If you decide to host jstree yourself - the files are located in the `dist` folder. You can safely ignore the `dist/libs` folder._

---

### Populating a tree using HTML

Now we are all set to create a tree, inline HTML is the easiest option (suitable for menus). All you need to do is select a node (using a jQuery selector) and invoke the `.jstree()` function to let jstree know you want to render a tree inside the selected node. `$.jstree.create(element)` can be used too.

```html
<div id="container">
  <ul>
    <li>Root node
      <ul>
        <li>Child node 1</li>
        <li>Child node 2</li>
      </ul>
    </li>
  </ul>
</div>
<script>
$(function() {
  $('#container').jstree();
});
</script>
```

[view result](http://jsfiddle.net/vakata/2kwkh2uL/)

_You can add a few options when rendering a node using a data-attribute (note the quotes):_
```html
<li data-jstree='{ "selected" : true, "opened" : true }'>Root node ...
```

---

### Populating a tree using an array (or JSON)

Building trees from HTML is easy, but it is not very flexible, inline JS data is a better option:

```html
<div id="container"></div>
<script>
$(function() {
  $('#container').jstree({
    'core' : {
      'data' : [
        { "text" : "Root node", "children" : [
            { "text" : "Child node 1" },
            { "text" : "Child node 2" }
          ]
        }
      ]
    }
  });
});
</script>
```

[view result](http://jsfiddle.net/vakata/2kwkh2uL/4478/)

Unlike the previous simple HTML example, this time the `.jstree()` function accepts a config object.

For now it is important to note that jstree will try to parse any data you specify in the  `core.data` key and use it to create a tree. As seen in the previous example, if this key is missing jstree will try to parse the inline HTML of the container.

#### The required JSON format

The data you use must be in a specific format, each branch of the tree is represented by an object, which must at least have a `text` key. The `children` key can be used to add children to the branch, it should be an array of objects.

_Keep in mind, you can use a simple string instead of an object if all you need is node with the given text, the above data can be written as:_

```js
[ { "text" : "Root node", "children" : [ "Child node 1", "Child node 2" ] } ]
```

There are other available options for each node, only set them if you need them like:

 * `id` - makes if possible to identify a node later (will also be used as a DOM ID of the `LI` node). _Make sure you do not repeat the same ID in a tree instance (that would defeat its purpose of being a unique identifier and may cause problems for jstree)_.
 * `icon` - a string which will be used for the node's icon - this can either be a path to a file, or a className (or list of classNames), which you can style in your CSS (font icons also work).
 * `data` - this can be anything you want - it is metadata you want attached to the node - you will be able to access and modify it any time later - it has no effect on the visuals of the node.
 * `state` - an object specifyng a few options about the node:
   - `selected` - if the node should be initially selected
   - `opened` - if the node should be initially opened
   - `disabled` - if the node should be disabled
   - `checked` - __checkbox plugin specific__ - if the node should be checked (only used when `tie_selection` is `false`, which you should only do if you really know what you are doing)
   - `undetermined` - __checkbox plugin specific__ - if the node should be rendered in undetermined state (only used with lazy loading and when the node is not yet loaded, otherwise this state is automatically calculated).
 * `type` - __types plugin specific__ - the type of the nodes (should be defined in the types config), if not set `"default"` is assumed.
 * `li_attr` - object of values which will be used to add HTML attributes on the resulting `LI` DOM node.
 * `a_attr` - object of values which will be used to add HTML attributes on the resulting `A` node.

Here is a new demo with some of those properties set:

```html
<div id="container"></div>
<script>
$(function() {
  $('#container').jstree({
    'core' : {
      'data' : [
          {
              "text" : "Root node",
              "state" : {"opened" : true },
              "children" : [
                  {
                    "text" : "Child node 1",
                    "state" : { "selected" : true },
                    "icon" : "glyphicon glyphicon-flash"
                  },
                  { "text" : "Child node 2", "state" : { "disabled" : true } }
              ]
        }
      ]
    }
  });
});
</script>
```

[view result](http://jsfiddle.net/vakata/2kwkh2uL/4479/)

---

### Populating the tree using AJAX

Building off of the previous example, let's see how to have jstree make AJAX requests for you.

```html
<div id="container"></div>
<script>
$(function() {
  $('#container').jstree({
    'core' : {
      'data' : {
        "url" : "//www.jstree.com/fiddle/",
        "dataType" : "json" // needed only if you do not supply JSON headers
      }
    }
  });
});
</script>
```

The server response is:
```json
[{
  "id":1,"text":"Root node","children":[
    {"id":2,"text":"Child node 1"},
    {"id":3,"text":"Child node 2"}
  ]
}]
```

[view result](http://jsfiddle.net/vakata/2kwkh2uL/4480/)

Instead of a JS array, you can set `core.data` to a [jQuery AJAX config](http://api.jquery.com/jQuery.ajax/). 
jsTree will hit that URL, and provided you return properly formatted JSON it will be displayed.

_If you cannot provide proper JSON headers, set `core.data.dataType` to `"json"`._

The ids in the server response make it possible to identify nodes later (which we will see in the next few demos), but they are not required.

__WHEN USING IDS MAKE SURE THEY ARE UNIQUE INSIDE A PARTICULAR TREE__

---

### Populating the tree using AJAX and lazy loading nodes

Lazy loading means nodes will be loaded when they are needed. Imagine you have a huge amount of nodes you want to show, but loading them with a single request is way too much traffic. Lazy loading makes it possible to load nodes on the fly - jstree will perform AJAX requests as the user browses the tree.

Here we take our previous example, and lazy load the "Child node 1" node.

```html
<div id="container"></div>
<script>
$(function() {
  $('#container').jstree({
    'core' : {
      'data' : {
        "url" : "//www.jstree.com/fiddle/?lazy",
        "data" : function (node) {
          return { "id" : node.id };
        }
      }
    }
  });
});
</script>
```

The initial server response is:
```json
[{
  "id":1,"text":"Root node","children":[
    {"id":2,"text":"Child node 1","children":true},
    {"id":3,"text":"Child node 2"}
  ]
}]
```

[view result](http://jsfiddle.net/vakata/2kwkh2uL/4481/)

Now to focus on what is different. First off the `"data"` config option of the data object. If you check with jQuery, it is supposed to be a string or an object. But jstree makes it possible to set a function.

Each time jstree needs to make an AJAX call this function will be called and will receive a single parameter - the node that is being loaded. The return value of this function will be used as the actual `"data"` of the AJAX call. To understand better open up the demo and see the requests go off in the console.

You will notice that the first request goes off to:
`http://www.jstree.com/fiddle?lazy&id=#`
`#` is the special ID that the function receives when jstree needs to load the root nodes.

Now go ahead and open the root node - two children will be shown, but no request will be made - that is because we loaded those children along with the first request.

Onto the next difference - "Child node 1" appears closed - that is because in the data we supplied `true` as the `"children"` property of this node (you can see it in the server response). This special value indicated to jstree, that it has to lazy load the "Child node 1" node.

Proceed and open this node - you will see a next request fire off to:
`http://www.jstree.com/fiddle?lazy&id=2`
ID is set to `2` because the node being loaded has an ID of `2`, and we have configured jstree to send the node ID along with the AJAX request (the `data` function).

The server response is:
```json
["Child node 3","Child node 4"]
```

_You can also set `"url"` to a function and it works exactly as with `"data"` - each time a request has to be made, jstree will invoke your function and the request will go off to whatever you return in this function. This is useful when dealing with URLs like: `http://example.com/get_children/1`._

### Populating the tree using a callback function

Sometimes you may not want jsTree to make AJAX calls for you - you might want to make them yourself, or use some other method of populating the tree. In that case you can use a callback function.

```html
<div id="container"></div>
<script>
$(function() {
  $('#container').jstree({
    'core' : {
      'data' : function (node, cb) {
        if(node.id === "#") {
          cb([{"text" : "Root", "id" : "1", "children" : true}]);
        }
        else {
          cb(["Child"]);
        }
      }
    }
  });
});
</script>
```

[view result](http://jsfiddle.net/vakata/2kwkh2uL/4482/)

As you can see your function will receive two arguments - the node whose children need to be loaded and a callback function to call with the data once you have it. The data follows the same familiar JSON format and lazy loading works just as with AJAX (as you can see in the above example).

---

## Working with events

jstree provides a lot of events to let you know something happened with the tree. The events are the same regardless of how you populate the tree.
Let's use the most basic event `changed` - it fires when selection on the tree changes:

```html
<div id="container"></div>
<script>
$(function() {
  $('#container').jstree({
    'core' : {
      'data' : [
        {"id" : 1, "text" : "Node 1"},
        {"id" : 2, "text" : "Node 2"},
      ]
    }
  });
  $('#container').on("changed.jstree", function (e, data) {
    console.log("The selected nodes are:");
    console.log(data.selected);
  });
});
</script>
```

[view result](http://jsfiddle.net/vakata/2kwkh2uL/4483/)

All jstree events fire in a special `".jstree"` namespace - this is why we listen for `"changed.jstree"`. The handler itself receives one additional parameter - it will be populated with all you need to know about the event that happened. In this case `data.selected` is an array of selected node IDs (please note, that if you have not specified IDs they will be autogenerated).

Let's extend this a bit and log out the text of the node instead of the ID.

```js
$('#container').on("changed.jstree", function (e, data) {
  console.log(data.instance.get_selected(true)[0].text);
  console.log(data.instance.get_node(data.selected[0]).text);
});
```

The two rows above achieve exactly the same thing - get the text of the first selected node.

In the `data` argument object you will always get an `instance` key - that is a reference to the tree instance, so that you can easily invoke methods.

__All available functions and events are documented in the API docs__

---

## Interacting with the tree using the API

We scratched the surface on interacting with the tree in the previous example. Let's move on to obtaining an instance and calling a method on this instance:

```html
<button>Select node 1</button>
<div id="container"></div>
<script>
$(function() {
  $('#container').jstree({
    'core' : {
      'data' : [
        {"id" : 1, "text" : "Node 1"},
        {"id" : 2, "text" : "Node 2"},
      ]
    }
  });
  $('button').on("click", function () {
    var instance = $('#container').jstree(true);
    instance.deselect_all();
    instance.select_node('1');
  });
});
</script>
```

[view result](http://jsfiddle.net/vakata/2kwkh2uL/4484/)

The above example shows how to obtain a reference to a jstree instance (again with a selector, but this time instead of a config, we pass a boolean `true`), and call a couple of methods - the latter one is selecting a node by its ID.

Methods can also be invoked like this:

```js
$('#container').jstree("select_node", "1");
```

__All available functions and events are documented in the API docs__

## More on configuration

We already covered the config object in general (when we specified inline & AJAX data sources).

```js
$("#tree").jstree({ /* config object goes here */ });
```

Each key in the config object corresponds to a plugin, and the value of that key is the configuration for that plugin. There are also two special keys `"core"` and `"plugins"`:
 * `"core"` stores the core configuration options
 * `"plugins"` is an array of plugin names (strings) you want active on the instance

When configuring you only need to set values that you want to be different from the defaults.

__All config options and defaults are documented in the API docs__

```js
$("#tree").jstree({
  "core" : { // core options go here
    "multiple" : false, // no multiselection
    "themes" : {
      "dots" : false // no connecting dots between dots
    }
  },
  "plugins" : ["state"] // activate the state plugin on this instance
});
```

[view result](http://jsfiddle.net/vakata/2kwkh2uL/4485/)

We will cover all plugins further down.

__Keep in mind by default all modifications to the structure are prevented - that means drag'n'drop, create, rename, delete will not work unless you enable them.__

```js
$("#tree").jstree({
  "core" : {
    "check_callback" : true, // enable all modifications
  },
  "plugins" : ["dnd","contextmenu"]
});
```

[view result](http://jsfiddle.net/vakata/2kwkh2uL/4486/)

`"core.check_callback"` can also be set to a function, that will be invoked every time a modification is about to happen (or when jstree needs to check if a modification is possible). If you return `true` the operation will be allowed, a value of `false` means it will not be allowed. The possible operation you can expect are `create_node`, `rename_node`, `delete_node`, `move_node` and `copy_node`. The `more` parameter will contain various information provided by the plugin that is invoking the check. For example the DND plugin will provide an object containing information about the move or copy operation that is being checked - is it a multi tree operation, which node is currently hovered, where the insert arrow is pointing - before, after or inside, etc.

```js
$("#tree").jstree({
  "core" : {
    "check_callback" : function (operation, node, parent, position, more) {
      if(operation === "copy_node" || operation === "move_node") {
        if(parent.id === "#") {
          return false; // prevent moving a child above or below the root
        }
      },
      return true; // allow everything else
    }
  },
  "plugins" : ["dnd","contextmenu"]
});
```

[view result](http://jsfiddle.net/vakata/2kwkh2uL/4487/)

The `more` parameter you receive contains other information related to the check being performed.

__For example__: `move_node` & `copy_node` checks will fire repeatedly while the user drags a node, if the check was triggered by the `dnd` plugin `more` will contain a `dnd` key, which will be set to `true`.
You can check for `more.dnd` and only perform a certain action if `dnd` triggered the check.
If you only want to perform an operation when a node is really about to be dropped check for `more.core`.

## Plugins

jsTree comes with a few plugin bundled, but they will only modify your tree if you activate them using the `"plugins"` config option. Here is a brief description of each plugin. You can read more on the available config options for each plugin in the API docs.

### checkbox
Renders a checkbox icon in front of each node, making multiselection easy. It also has a "tri-state" option, meaning a node with some of its children checked will get a "square" icon.

_Keep in mind that if any sort of cascade is enabled, disabled nodes may be checked too (not by themselves, but for example when a parent of a disabled node is checked and selection is configured to cascade down)._

```js
$("#tree").jstree({
  "plugins" : ["checkbox"]
});
```

[view result](http://jsfiddle.net/vakata/2kwkh2uL/4488/)

### contextmenu
Makes it possible to right click nodes and shows a list of configurable actions in a menu.

```js
$("#tree").jstree({
  "core" : { "check_callback" : true }, // so that modifying operations work
  "plugins" : ["contextmenu"]
});
```

[view result](http://jsfiddle.net/vakata/2kwkh2uL/4489/)

### dnd
Makes it possible to drag and drop tree nodes and rearrange the tree.

```js
$("#tree").jstree({
  "core" : { "check_callback" : true }, // so that operations work
  "plugins" : ["dnd"]
});
```

[view result](http://jsfiddle.net/vakata/2kwkh2uL/4490/)

### massload
Makes it possible to load multiple nodes in a single go (for a lazy loaded tree).

```js
$("#tree").jstree({
  "core" : {
    "data" : { .. AJAX config .. }
  },
  "massload" : {
    "url" : "/some/path",
    "data" : function (nodes) {
      return { "ids" : nodes.join(",") };
    }
  },
  "plugins" : [ "massload", "state" ]
});
```

### search
Adds the possibility to search for items in the tree and show only matching nodes. It also has AJAX / callback hooks, so that search will work on lazy loaded trees too.

```html
<form id="s">
  <input type="search" id="q" />
  <button type="submit">Search</button>
</form>
<script>
$("#container").jstree({
  "plugins" : ["search"]
});
$("#s").submit(function(e) {
  e.preventDefault();
  $("#container").jstree(true).search($("#q").val());
});
</script>
```

[view result](http://jsfiddle.net/vakata/2kwkh2uL/4491/)

### sort
Automatically arranges all sibling nodes according to a comparison config option function, which defaults to alphabetical order.

```js
$("#tree").jstree({
  "plugins" : ["sort"]
});
```

[view result](http://jsfiddle.net/vakata/2kwkh2uL/4492/)

### state
Saves all opened and selected nodes in the user's browser, so when returning to the same tree the previous state will be restored.

```js
$("#tree").jstree({
  // the key is important if you have multiple trees in the same domain
  "state" : { "key" : "state_demo" },
  "plugins" : ["state"]
});
```

[view result](http://jsfiddle.net/vakata/2kwkh2uL/4493/)

### types
Makes it possible to add a "type" for a node, which means to easily control nesting rules and icon for groups of nodes instead of individually. To set a node type add a type property to the node structure.

```js
$("#tree").jstree({
  "types" : {
    "default" : {
      "icon" : "glyphicon glyphicon-flash"
    },
    "demo" : {
      "icon" : "glyphicon glyphicon-ok"
    }
  },
  "plugins" : ["types"]
});
```

[view result](http://jsfiddle.net/vakata/2kwkh2uL/4494/)

### unique
Enforces that no nodes with the same name can coexist as siblings - prevents renaming and moving nodes to a parent, which already contains a node with the same name.

```js
$("#tree").jstree({
  "plugins" : ["unique"]
});
```

[view result](http://jsfiddle.net/vakata/2kwkh2uL/4495/)

### wholerow
Makes each node appear block level which makes selection easier. May cause slow down for large trees in old browsers.

```js
$("#tree").jstree({
  "plugins" : ["wholerow"]
});
```

[view result](http://jsfiddle.net/vakata/2kwkh2uL/4496/)

### More plugins
If you create your own plugin (or download a 3rd party one) you must include its source on the page and list its name in the `"plugins"` config array.

```js
// conditional select
(function ($, undefined) {
  "use strict";
  $.jstree.defaults.conditionalselect = function () { return true; };
  $.jstree.plugins.conditionalselect = function (options, parent) {
    this.activate_node = function (obj, e) {
      if(this.settings.conditionalselect.call(this, this.get_node(obj))) {
        parent.activate_node.call(this, obj, e);
      }
    };
  };
})(jQuery);
$("#tree").jstree({
  "conditionalselect" : function (node) {
    return node.text === "Root node" ? false : true;
  },
  "plugins" : ["conditionalselect"]
});
```

[view result](http://jsfiddle.net/vakata/2kwkh2uL/4497/)

As seen here when creating a plugin you can define a default config, add your own functions to jstree, or override existing ones while maintaining the ability to call the overridden function.

## PHP demos moved to new repository
https://github.com/vakata/jstree-php-demos

## License & Contributing

_Please do NOT edit files in the "dist" subdirectory as they are generated via grunt. You'll find source code in the "src" subdirectory!_

If you want to you can always [donate a small amount][paypal] to help the development of jstree.

[paypal]: https://www.paypal.com/cgi-bin/webscr?cmd=_xclick&business=paypal@vakata.com&currency_code=USD&amount=&return=http://jstree.com/donation&item_name=Buy+me+a+coffee+for+jsTree

Copyright (c) 2014 Ivan Bozhanov (http://vakata.com)

Licensed under the [MIT license](http://www.opensource.org/licenses/mit-license.php).
# SQLite compiled to javascript
[![Build Status](https://travis-ci.org/kripken/sql.js.svg?branch=master)](http://travis-ci.org/kripken/sql.js) [![CDNJS version](https://img.shields.io/cdnjs/v/sql.js.svg)](https://cdnjs.com/libraries/sql.js)

For the impatients, try the demo here: http://kripken.github.io/sql.js/examples/GUI

*sql.js* is a port of [SQLite](http://sqlite.org/about.html) to Webassembly, by compiling the SQLite C code with [Emscripten](http://kripken.github.io/emscripten-site/docs/introducing_emscripten/about_emscripten.html). It uses a [virtual database file stored in memory](https://kripken.github.io/emscripten-site/docs/porting/files/file_systems_overview.html), and thus **doesn't persist the changes** made to the database. However, it allows you to **import** any existing sqlite file, and to **export** the created database as a [javascript typed array](https://developer.mozilla.org/en-US/docs/Web/JavaScript/Typed_arrays).

There are no C bindings or node-gyp compilation here, sql.js is a simple javascript file, that can be used like any traditional javascript library. If you are building a native application in javascript (using Electron for instance), or are working in node.js, you will likely prefer to use [a native binding of SQLite to javascript](https://www.npmjs.com/package/sqlite3).

SQLite is public domain, sql.js is MIT licensed.

Sql.js predates WebAssembly, and thus started as an [asm.js](https://en.wikipedia.org/wiki/Asm.js) project. It still supports asm.js for backwards compatability.

## Version of binaries
Sql.js was last built with:
Emscripten version 1.38.30 (2019-04-16) [Release History](https://emscripten.org/docs/introducing_emscripten/release_notes.html)
SqlLite version: 3.28.0 (2019-04-16) [Release History](https://www.sqlite.org/changes.html)

## Documentation
A [full documentation](http://kripken.github.io/sql.js/documentation/#http://kripken.github.io/sql.js/documentation/class/Database.html) generated from comments inside the source code, is available.

## Usage

```javascript
var initSqlJs = require('sql-wasm.js');
// or if you are in a browser:
//var initSqlJs = window.initSqlJs;

initSqlJs().then(function(SQL){

  // Create a database
  var db = new SQL.Database();
  // NOTE: You can also use new SQL.Database(data) where
  // data is an Uint8Array representing an SQLite database file

  // Execute some sql
  sqlstr = "CREATE TABLE hello (a int, b char);";
  sqlstr += "INSERT INTO hello VALUES (0, 'hello');"
  sqlstr += "INSERT INTO hello VALUES (1, 'world');"
  db.run(sqlstr); // Run the query without returning anything

  var res = db.exec("SELECT * FROM hello");
  /*
  [
    {columns:['a','b'], values:[[0,'hello'],[1,'world']]}
  ]
  */

  // Prepare an sql statement
  var stmt = db.prepare("SELECT * FROM hello WHERE a=:aval AND b=:bval");

  // Bind values to the parameters and fetch the results of the query
  var result = stmt.getAsObject({':aval' : 1, ':bval' : 'world'});
  console.log(result); // Will print {a:1, b:'world'}

  // Bind other values
  stmt.bind([0, 'hello']);
  while (stmt.step()) console.log(stmt.get()); // Will print [0, 'hello']

  // You can also use javascript functions inside your SQL code
  // Create the js function you need
  function add(a, b) {return a+b;}
  // Specifies the SQL function's name, the number of it's arguments, and the js function to use
  db.create_function("add_js", add);
  // Run a query in which the function is used
  db.run("INSERT INTO hello VALUES (add_js(7, 3), add_js('Hello ', 'world'));"); // Inserts 10 and 'Hello world'

  // free the memory used by the statement
  stmt.free();
  // You can not use your statement anymore once it has been freed.
  // But not freeing your statements causes memory leaks. You don't want that.

  // Export the database to an Uint8Array containing the SQLite database file
  var binaryArray = db.export();
});

```

## Demo
There are a few examples [available here](https://kripken.github.io/sql.js/index.html). The most full-featured is the [Sqlite Interpreter](https://kripken.github.io/sql.js/examples/GUI/index.html).

## Examples
The test files provide up to date example of the use of the api.
### Inside the browser
#### Example **HTML** file:
```html
<meta charset="utf8" />
<html>
  <script src='/dist/sql-wasm.js'></script>
  <script>
    config = {
      locateFile: url => `/dist/${filename}` 
    }
    // The `initSqlJs` function is globally provided by all of the main dist files if loaded in the browser.
    // We must specify this locateFile function if we are loading a wasm file from anywhere other than the current html page's folder.
    initSqlJs(config).then(function(SQL){
      //Create the database
      var db = new SQL.Database();
      // Run a query without reading the results
      db.run("CREATE TABLE test (col1, col2);");
      // Insert two rows: (1,111) and (2,222)
      db.run("INSERT INTO test VALUES (?,?), (?,?)", [1,111,2,222]);
  
      // Prepare a statement
      var stmt = db.prepare("SELECT * FROM test WHERE col1 BETWEEN $start AND $end");
      stmt.getAsObject({$start:1, $end:1}); // {col1:1, col2:111}
  
      // Bind new values
      stmt.bind({$start:1, $end:2});
      while(stmt.step()) { //
        var row = stmt.getAsObject();
        console.log('Here is a row: ' + JSON.stringify(row));
      }
    });
  </script>
  <body>
    Output is in Javscript console
  </body>
</html>
```

#### Creating a database from a file choosen by the user
`SQL.Database` constructor takes an array of integer representing a database file as an optional parameter.
The following code uses an HTML input as the source for loading a database:
```javascript
dbFileElm.onchange = () => {
  var f = dbFileElm.files[0];
  var r = new FileReader();
  r.onload = function() {
    var Uints = new Uint8Array(r.result);
    db = new SQL.Database(Uints);
  }
  r.readAsArrayBuffer(f);
}
```
See : http://kripken.github.io/sql.js/examples/GUI/gui.js

#### Loading a database from a server

```javascript
var xhr = new XMLHttpRequest();
// For example: https://github.com/lerocha/chinook-database/raw/master/ChinookDatabase/DataSources/Chinook_Sqlite.sqlite
xhr.open('GET', '/path/to/database.sqlite', true);
xhr.responseType = 'arraybuffer';

xhr.onload = e => {
  var uInt8Array = new Uint8Array(this.response);
  var db = new SQL.Database(uInt8Array);
  var contents = db.exec("SELECT * FROM my_table");
  // contents is now [{columns:['col1','col2',...], values:[[first row], [second row], ...]}]
};
xhr.send();
```
See: https://github.com/kripken/sql.js/wiki/Load-a-database-from-the-server


### Use from node.js

`sql.js` is [hosted on npm](https://www.npmjs.org/package/sql.js). To install it, you can simply run `npm install sql.js`.
Alternatively, you can simply download `sql-wasm.js` and `sql-wasm.wasm`, from the download link below.

#### read a database from the disk:
```javascript
var fs = require('fs');
var initSqlJs = require('sql-wasm.js');
var filebuffer = fs.readFileSync('test.sqlite');
 
initSqlJs().then(function(SQL){
  // Load the db
  var db = new SQL.Database(filebuffer);
});

```

#### write a database to the disk
You need to convert the result of `db.export` to a buffer
```javascript
var fs = require("fs");
// [...] (create the database)
var data = db.export();
var buffer = new Buffer(data);
fs.writeFileSync("filename.sqlite", buffer);
```

See : https://github.com/kripken/sql.js/blob/master/test/test_node_file.js

### Use as web worker
If you don't want to run CPU-intensive SQL queries in your main application thread,
you can use the *more limited* WebWorker API.

You will need to download [dist/worker.sql-wasm.js](dist/worker.sql-wasm.js) [dist/worker.sql-wasm.wasm](dist/worker.sql-wasm.wasm).

Example:
```html
<script>
  var worker = new Worker("/dist/worker.sql-wasm.js");
  worker.onmessage = () => {
    console.log("Database opened");
    worker.onmessage = event => {
      console.log(event.data); // The result of the query
    };
	
    worker.postMessage({
      id: 2,
      action: 'exec',
      sql: 'SELECT * FROM test'
    });
  };

  worker.onerror = e => console.log("Worker error: ", e);
  worker.postMessage({
    id:1,
    action:'open',
    buffer:buf, /*Optional. An ArrayBuffer representing an SQLite Database file*/
  });
</script>
```

See [examples/GUI/gui.js](examples/GUI/gui.js) for a full working example.

## Flavors/versions Targets/Downloads

This library includes both WebAssembly and asm.js versions of Sqlite. (WebAssembly is the newer, preferred way to compile to Javascript, and has superceded asm.js. It produces smaller, faster code.) Asm.js versions are included for compatibility.

## Upgrading from 0.x to 1.x

Version 1.0 of sql.js must be loaded asynchronously, whereas asm.js was able to be loaded synchronously. 

So in the past, you would:
```html
<script src='js/sql.js'></script>
<script>
  var db = new SQL.Database();
  //...
</script>
```
or:
```javascript
var SQL = require('sql.js');
var db = new QL.Database();
//...
```

Version 1.x:
```html
<script src='dist/sql-wasm.js'></script>
<script>
  initSqlJs({ locateFile: filename => `/dist/${filename}` }).then(function(SQL){
    var db = new SQL.Database();
    //...
  });
</script>
```
or:
```javascript
var initSqlJs = require('sql-wasm.js');
initSqlJs().then(function(SQL){
  var db = new SQL.Database();
  //...
});
```

`NOTHING` is now a reserved word in SQLite, whereas previously it was not. This could cause errors like `Error: near "nothing": syntax error`

### Downloading/Using: ###
Although asm.js files were distributed as a single Javascript file, WebAssembly libraries are most efficiently distributed as a pair of files, the `.js`  loader and the `.wasm` file, like [dist/sql-wasm.js]([dist/sql-wasm.js]) and [dist/sql-wasm.wasm]([dist/sql-wasm.wasm]). The `.js` file is reponsible for wrapping/loading the `.wasm` file. 




## Versions of sql.js included in `dist/`
 - `sql-wasm.js` : The Web Assembly version of Sql.js. Minified and suitable for production. Use this. If you use this, you will need to include/ship `sql-wasm.wasm` as well.
 - `sql-wasm-debug.js` : The Web Assembly, Debug version of Sql.js. Larger, with assertions turned on. Useful for local development. You will need to include/ship `sql-wasm-debug.wasm` if you use this.
 - `sql-asm.js` : The older asm.js version of Sql.js. Slower and larger. Provided for compatiblity reasons.
 - `sql-asm-memory-growth.js` : Asm.js doesn't allow for memory to grow by default, because it is slower and de-optimizes. If you are using sql-asm.js and you see this error (`Cannot enlarge memory arrays`), use this file.
 - `sql-asm-debug.js` : The _Debug_ asm.js version of Sql.js. Use this for local development.
 - `worker.*` - Web Worker versions of the above libraries. More limited API. See [examples/GUI/gui.js](examples/GUI/gui.js) for a good example of this.

## Compiling

- Install the EMSDK, [as described here](https://kripken.github.io/emscripten-site/docs/getting_started/downloads.html)
- Run `npm run rebuild`

# GeoPackage JS

GeoPackage JS is an implementation of the OGC GeoPackage spec.  This library works in both the browser and Node 4+.

### Demo ###
[GeoPackage JS Demo Page](http://ngageoint.github.io/geopackage-js/)

Cloning this repository and opening the docs/index.html in your browser will run the demo locally.

### Installation ###

[![Build Status](https://travis-ci.org/ngageoint/geopackage-js.svg?branch=master)](https://travis-ci.org/ngageoint/geopackage-js)
[![NPM](https://img.shields.io/npm/v/@ngageoint/geopackage.svg)](https://www.npmjs.com/package/@ngageoint/geopackage)
[![Coverage Status](https://coveralls.io/repos/github/ngageoint/geopackage-js/badge.svg)](https://coveralls.io/github/ngageoint/geopackage-js)

```sh
$ npm install @ngageoint/geopackage
```

#### GeoPackage JS Library ####

The [GeoPackage Libraries](http://ngageoint.github.io/GeoPackage/) were developed at the [National Geospatial-Intelligence Agency (NGA)](http://www.nga.mil/) in collaboration with [BIT Systems](http://www.bit-sys.com/). The government has "unlimited rights" and is releasing this software to increase the impact of government investments by providing developers with the opportunity to take things in new directions. The software use, modification, and distribution rights are stipulated within the [MIT license](http://choosealicense.com/licenses/mit/).

### Pull Requests ###
If you'd like to contribute to this project, please make a pull request. We'll review the pull request and discuss the changes. All pull request contributions to this project will be released under the MIT license.

Software source code previously released under an open source license and then modified by NGA staff is considered a "joint work" (see 17 USC § 101); it is partially copyrighted, partially public domain, and as a whole is protected by the copyrights of the non-government authors and must be released according to the terms of the original open source license.

### About ###

[GeoPackage JS](https://github.com/ngageoint/geopackage-js) is a [GeoPackage Library](http://ngageoint.github.io/GeoPackage/) JavaScript implementation of the Open Geospatial Consortium [GeoPackage](http://www.geopackage.org/) [spec](http://www.geopackage.org/spec/).  It is listed as an [OGC GeoPackage Implementation](http://www.geopackage.org/#implementations_nga) by the National Geospatial-Intelligence Agency.

The GeoPackage JavaScript library currently provides the ability to read GeoPackage files.  This library works both in the browser and in Node.  In the browser tiles are rendered using HTML5 Canvas and GeoPackages are read using [sql.js](https://github.com/kripken/sql.js/).  In Node tiles are rendered  [PureImage](https://github.com/joshmarinacci/node-pureimage) and GeoPackages are read using [node-sqlite3](https://github.com/mapbox/node-sqlite3).

### Changelog

##### 2.1.0

- Implementation of the Feature Style Extension and Contents ID Extension

##### 2.0.8

- Checks for Electron when returning a tile creator

##### 2.0

- All new API utilizing Promises

##### 1.1.4

- Adds a method to retrieve tiles in EPSG:4326

##### 1.1.3

- Fixes issue #115

##### 1.1.2

- fix case where GeoPackage Zoom does not correspond to the web map zoom

##### 1.1.1

- fix more instances of proj4 bug for react
- fixed tile generation for images with different x and y pixel densities

##### 1.1.0

- accept pull request adding support for react
- fix bug with projected tiles that spanned the date line

##### 1.0.25

- ensure we use proj4 2.4.3 instead of 2.4.4

##### 1.0.22

- Fixed bug where querying for indexed features only returned the geometry instead of the entire feature

##### 1.0.19

- Remove dependency on Lwip

### Usage ###

View examples using [Bower](https://github.com/ngageoint/geopackage-js/tree/master/docs/bower) and [Browserify](https://github.com/ngageoint/geopackage-js/tree/master/docs)

View the latest [docs](http://ngageoint.github.io/geopackage-js/jsdoc/module-geoPackage-GeoPackage.html) (currently being updated).

#### Browser Usage ####
```javascript

// attach this method to a file input onchange event
window.loadGeoPackage = function(files) {
  var f = files[0];
  var r = new FileReader();
  r.onload = function() {
    var array = new Uint8Array(r.result);
    loadByteArray(array);
  }
  r.readAsArrayBuffer(f);
}

function loadByteArray(array, callback) {
  var db = new SQL.Database(array);
  GeoPackageConnection.connectWithDatabase(db, function(err, connection) {
    var geoPackage = new GeoPackage('', '', connection);

    // Now you can operate on the GeoPackage

    // get the tile table names
    geoPackage.getTileTables(function(err, tileTableNames) {
      // tileTableNames is an array of all tile table names

      // get the info for the first table
      geoPackage.getTileDaoWithTableName(tileTableNames[0], function(err, tileDao) {
        geoPackage.getInfoForTable(tileDao, function(err, info) {
          // do something with the tile table info
        });

        // draw a tile into a canvas for an XYZ tile
        var canvas = canvasFromSomewhere;
        var gpr = new GeoPackageTileRetriever(tileDao, 256, 256);
        var x = 0;
        var y = 0;
        var zoom = 0;

        console.time('Draw tile ' + x + ', ' + y + ' zoom: ' + zoom);
        gpr.drawTileIn(x, y, zoom, canvas, function() {
          console.timeEnd('Draw tile ' + x + ', ' + y + ' zoom: ' + zoom);
        });

        // or get a tile base64 data URL for an XYZ tile
        gpr.getTile(x, y, zoom, function(err, tileBase64DataURL) {
          console.log('got the base64 data url');
        });

        // or get a tile from a GeoPackage tile column and tile row
        tileDao.queryForTile(tileColumn, tileRow, zoom, function(err, tile) {
          var tileData = tile.getTileData();  // the raw bytes from the GeoPackage
        });

      });
    });

    // get the feature table names
    geoPackage.getFeatureTables(function(err, featureTableNames) {
      // featureTableNames is an array of all feature table names

      // get the info for the first table
      geoPackage.getFeatureDaoWithTableName(featureTableNames[0], function(err, featureDao) {
        geoPackage.getInfoForTable(featureDao, function(err, info) {
          // do something with the feature table info
        });

        // query for all features
        featureDao.queryForEach(function(err, row, rowDone) {
          var feature = featureDao.getFeatureRow(row);
          var geometry = currentRow.getGeometry();
          if (geometry) {
            var geom = geometry.geometry;
            var geoJson = geometry.geometry.toGeoJSON();

            geoJson.properties = {};
            for (var key in feature.values) {
              if(feature.values.hasOwnProperty(key) && key != feature.getGeometryColumn().name) {
                var column = info.columnMap[key];
                geoJson.properties[column.displayName] = currentRow.values[key];
              }
            }
          }
          rowDone();
        });
      });
    });
  });
}

```

#### NodeJS Usage ####

```javascript
var GeoPackageAPI = require('@ngageoint/geopackage')
  , GeoPackageManager = GeoPackageAPI.GeoPackageManager
  , GeoPackageConnection = GeoPackageAPI.GeoPackageConnection
  , GeoPackageTileRetriever = GeoPackageAPI.GeoPackageTileRetriever;

GeoPackageAPI.open(filename, function(err, geoPackage) {

  // Now you can operate on the GeoPackage

  // get the tile table names
  geoPackage.getTileTables(function(err, tileTableNames) {
    // tileTableNames is an array of all tile table names

    // get the info for the first table
    geoPackage.getTileDaoWithTableName(tileTableNames[0], function(err, tileDao) {
      geoPackage.getInfoForTable(tileDao, function(err, info) {
        // do something with the tile table info
      });

      // draw a tile into a canvas for an XYZ tile
      var canvas = canvasFromSomewhere;
      var gpr = new GeoPackageTileRetriever(tileDao, 256, 256);
      var x = 0;
      var y = 0;
      var zoom = 0;

      console.time('Draw tile ' + x + ', ' + y + ' zoom: ' + zoom);
      gpr.drawTileIn(x, y, zoom, canvas, function() {
        console.timeEnd('Draw tile ' + x + ', ' + y + ' zoom: ' + zoom);
      });

      // or get a tile base64 data URL for an XYZ tile
      gpr.getTile(x, y, zoom, function(err, tileBase64DataURL) {
        console.log('got the base64 data url');
      });

      // or get a tile from a GeoPackage tile column and tile row
      tileDao.queryForTile(tileColumn, tileRow, zoom, function(err, tile) {
        var tileData = tile.getTileData();  // the raw bytes from the GeoPackage
      });

    });
  });

  // get the feature table names
  geoPackage.getFeatureTables(function(err, featureTableNames) {
    // featureTableNames is an array of all feature table names

    // get the info for the first table
    geoPackage.getFeatureDaoWithTableName(featureTableNames[0], function(err, featureDao) {
      geoPackage.getInfoForTable(featureDao, function(err, info) {
        // do something with the feature table info
      });

      // query for all features
      featureDao.queryForEach(function(err, row, rowDone) {
        var feature = featureDao.getFeatureRow(row);
        var geometry = currentRow.getGeometry();
        if (geometry) {
          var geom = geometry.geometry;
          var geoJson = geometry.geometry.toGeoJSON();

          geoJson.properties = {};
          for (var key in feature.values) {
            if(feature.values.hasOwnProperty(key) && key != feature.getGeometryColumn().name) {
              var column = info.columnMap[key];
              geoJson.properties[column.displayName] = currentRow.values[key];
            }
          }
        }
        rowDone();
      });
    });
  });
});

```
# MeshLine
Mesh replacement for ```THREE.Line```

Instead of using GL_LINE, it uses a strip of triangles billboarded. Some examples:

[![Demo](screenshots/demo.jpg)](https://www.clicktorelease.com/code/THREE.MeshLine/demo/index.html)
[![Graph](screenshots/graph.jpg)](https://www.clicktorelease.com/code/THREE.MeshLine/demo/graph.html)
[![Spinner](screenshots/spinner.jpg)](https://www.clicktorelease.com/code/THREE.MeshLine/demo/spinner.html)
[![SVG](screenshots/svg.jpg)](https://www.clicktorelease.com/code/THREE.MeshLine/demo/svg.html)
[![Shape](screenshots/shape.jpg)](https://www.clicktorelease.com/code/THREE.MeshLine/demo/shape.html)
[![Shape](screenshots/birds.jpg)](https://www.clicktorelease.com/code/THREE.MeshLine/demo/birds.html)

* [Demo](https://www.clicktorelease.com/code/THREE.MeshLine/demo/index.html): play with the different settings of materials
* [Graph](https://www.clicktorelease.com/code/THREE.MeshLine/demo/graph.html): example of using ```MeshLine``` to plot graphs
* [Spinner](https://www.clicktorelease.com/code/THREE.MeshLine/demo/spinner.html): example of dynamic ```MeshLine``` with texture
* [SVG](https://www.clicktorelease.com/code/THREE.MeshLine/demo/svg.html): example of ```MeshLine``` rendering SVG Paths
* [Shape](https://www.clicktorelease.com/code/THREE.MeshLine/demo/shape.html): example of ```MeshLine``` created from a mesh
* [Birds](https://www.clicktorelease.com/code/THREE.MeshLine/demo/birds.html): example of ```MeshLine.advance()``` by @caramelcode (Jared Sprague) and @mwcz (Michael Clayton)

### How to use ####

* Include script
* Create and populate a geometry
* Create a MeshLine and assign the geometry
* Create a MeshLineMaterial
* Use MeshLine and MeshLineMaterial to create a THREE.Mesh

#### Include the script

Include script after THREE is included
```js
<script src="THREE.MeshLine.js"></script>
```
or use npm to install it
```
npm i three.meshline
```
and include it in your code (don't forget to require three.js)
```js
var THREE = require( 'three' );
var MeshLine = require( 'three.meshline' );
```

##### Create and populate a geometry #####

First, create the list of vertices that will define the line. ```MeshLine``` accepts ```THREE.Geometry``` (looking up the ```.vertices``` in it) and ```Array```/```Float32Array```. ```THREE.BufferGeometry``` coming soon, and may be others like ```Array``` of ```THREE.Vector3```.

```js
var geometry = new THREE.Geometry();
for( var j = 0; j < Math.PI; j += 2 * Math.PI / 100 ) {
	var v = new THREE.Vector3( Math.cos( j ), Math.sin( j ), 0 );
	geometry.vertices.push( v );
}
```

##### Create a MeshLine and assign the geometry #####

Once you have that, you can create a new ```MeshLine```, and call ```.setGeometry()``` passing the vertices.

```js
var line = new MeshLine();
line.setGeometry( geometry );
```

Note: ```.setGeometry``` accepts a second parameter, which is a function to define the width in each point along the line. By default that value is 1, making the line width 1 * lineWidth.

```js
line.setGeometry( geometry, function( p ) { return 2; } ); // makes width 2 * lineWidth
line.setGeometry( geometry, function( p ) { return 1 - p; } ); // makes width taper
line.setGeometry( geometry, function( p ) { return 2 + Math.sin( 50 * p ); } ); // makes width sinusoidal
```

##### Create a MeshLineMaterial #####

A ```MeshLine``` needs a ```MeshLineMaterial```:

```js
var material = new MeshLineMaterial();
```

By default it's a white material of width 1 unit.

```MeshLineMaterial``` has several attributes to control the appereance of the ```MeshLine```:

* ```map``` - a ```THREE.Texture``` to paint along the line (requires ```useMap``` set to true)
* ```useMap``` - tells the material to use ```map``` (0 - solid color, 1 use texture)
* ```alphaMap``` - a ```THREE.Texture``` to use as alpha along the line (requires ```useAlphaMap``` set to true)
* ```useAlphaMap``` - tells the material to use ```alphaMap``` (0 - no alpha, 1 modulate alpha)
* ```repeat``` - THREE.Vector2 to define the texture tiling (applies to map and alphaMap - MIGHT CHANGE IN THE FUTURE)
* ```color``` - ```THREE.Color``` to paint the line width, or tint the texture with
* ```opacity``` - alpha value from 0 to 1 (requires ```transparent``` set to ```true```)
* ```alphaTest``` - cutoff value from 0 to 1
* ```dashArray``` - THREE.Vector2 to define the dashing (NOT IMPLEMENTED YET)
* ```resolution``` - ```THREE.Vector2``` specifying the canvas size (REQUIRED)
* ```sizeAttenuation``` - makes the line width constant regardless distance (1 unit is 1px on screen) (0 - attenuate, 1 - don't attenuate)
* ```lineWidth``` - float defining width (if ```sizeAttenuation``` is true, it's world units; else is screen pixels)
* ```near``` - camera near clip plane distance  (REQUIRED if ```sizeAttenuation``` set to false)
* ```far``` - camera far clip plane distance  (REQUIRED if ```sizeAttenuation``` set to false)

If you're rendering transparent lines or using a texture with alpha map, you should set ```depthTest``` to ```false```, ```transparent``` to ```true``` and ```blending``` to an appropriate blending mode, or use ```alphaTest```.

##### Use MeshLine and MeshLineMaterial to create a THREE.Mesh #####

Finally, we create a mesh and add it to the scene:

```js
var mesh = new THREE.Mesh( line.geometry, material ); // this syntax could definitely be improved!
scene.add( mesh );
```

### TODO ###

* Better miters
* Proper sizes
* Support for dashArray

### Support ###

Tested successfully on

* Chrome OSX, Windows, Android
* Firefox OSX, Windows, Anroid
* Safari OSX, iOS
* Internet Explorer 11 (SVG and Shape demo won't work because they use Promises)
* Opera OSX, Windows

### References ###

* [Drawing lines is hard](http://mattdesl.svbtle.com/drawing-lines-is-hard)
* [WebGL rendering of solid trails](http://codeflow.org/entries/2012/aug/05/webgl-rendering-of-solid-trails/)
* [Drawing Antialiased Lines with OpenGL](https://www.mapbox.com/blog/drawing-antialiased-lines/)

#### License ####

MIT licensed

Copyright (C) 2015-2016 Jaume Sanchez Elias, http://www.clicktorelease.com
# JSON5 – JSON for Humans

[![Build Status](https://travis-ci.org/json5/json5.svg)][Build Status]
[![Coverage
Status](https://coveralls.io/repos/github/json5/json5/badge.svg)][Coverage
Status]

The JSON5 Data Interchange Format (JSON5) is a superset of [JSON] that aims to
alleviate some of the limitations of JSON by expanding its syntax to include
some productions from [ECMAScript 5.1].

This JavaScript library is the official reference implementation for JSON5
parsing and serialization libraries.

[Build Status]: https://travis-ci.org/json5/json5

[Coverage Status]: https://coveralls.io/github/json5/json5

[JSON]: https://tools.ietf.org/html/rfc7159

[ECMAScript 5.1]: https://www.ecma-international.org/ecma-262/5.1/

## Summary of Features
The following ECMAScript 5.1 features, which are not supported in JSON, have
been extended to JSON5.

### Objects
- Object keys may be an ECMAScript 5.1 _[IdentifierName]_.
- Objects may have a single trailing comma.

### Arrays
- Arrays may have a single trailing comma.

### Strings
- Strings may be single quoted.
- Strings may span multiple lines by escaping new line characters.
- Strings may include character escapes.

### Numbers
- Numbers may be hexadecimal.
- Numbers may have a leading or trailing decimal point.
- Numbers may be [IEEE 754] positive infinity, negative infinity, and NaN.
- Numbers may begin with an explicit plus sign.

### Comments
- Single and multi-line comments are allowed.

### White Space
- Additional white space characters are allowed.

[IdentifierName]: https://www.ecma-international.org/ecma-262/5.1/#sec-7.6

[IEEE 754]: http://ieeexplore.ieee.org/servlet/opac?punumber=4610933

## Short Example
```js
{
  // comments
  unquoted: 'and you can quote me on that',
  singleQuotes: 'I can use "double quotes" here',
  lineBreaks: "Look, Mom! \
No \\n's!",
  hexadecimal: 0xdecaf,
  leadingDecimalPoint: .8675309, andTrailing: 8675309.,
  positiveSign: +1,
  trailingComma: 'in objects', andIn: ['arrays',],
  "backwardsCompatible": "with JSON",
}
```

## Specification
For a detailed explanation of the JSON5 format, please read the [official
specification](https://json5.github.io/json5-spec/).

## Installation
### Node.js
```sh
npm install json5
```

```js
const JSON5 = require('json5')
```

### Browsers
```html
<script src="https://unpkg.com/json5@^2.0.0/dist/index.min.js"></script>
```

This will create a global `JSON5` variable.

## API
The JSON5 API is compatible with the [JSON API].

[JSON API]:
https://developer.mozilla.org/en-US/docs/Web/JavaScript/Reference/Global_Objects/JSON

### JSON5.parse()
Parses a JSON5 string, constructing the JavaScript value or object described by
the string. An optional reviver function can be provided to perform a
transformation on the resulting object before it is returned.

#### Syntax
    JSON5.parse(text[, reviver])

#### Parameters
- `text`: The string to parse as JSON5.
- `reviver`: If a function, this prescribes how the value originally produced by
  parsing is transformed, before being returned.

#### Return value
The object corresponding to the given JSON5 text.

### JSON5.stringify()
Converts a JavaScript value to a JSON5 string, optionally replacing values if a
replacer function is specified, or optionally including only the specified
properties if a replacer array is specified.

#### Syntax
    JSON5.stringify(value[, replacer[, space]])
    JSON5.stringify(value[, options])

#### Parameters
- `value`: The value to convert to a JSON5 string.
- `replacer`: A function that alters the behavior of the stringification
  process, or an array of String and Number objects that serve as a whitelist
  for selecting/filtering the properties of the value object to be included in
  the JSON5 string. If this value is null or not provided, all properties of the
  object are included in the resulting JSON5 string.
- `space`: A String or Number object that's used to insert white space into the
  output JSON5 string for readability purposes. If this is a Number, it
  indicates the number of space characters to use as white space; this number is
  capped at 10 (if it is greater, the value is just 10). Values less than 1
  indicate that no space should be used. If this is a String, the string (or the
  first 10 characters of the string, if it's longer than that) is used as white
  space. If this parameter is not provided (or is null), no white space is used.
  If white space is used, trailing commas will be used in objects and arrays.
- `options`: An object with the following properties:
  - `replacer`: Same as the `replacer` parameter.
  - `space`: Same as the `space` parameter.
  - `quote`: A String representing the quote character to use when serializing
    strings.

#### Return value
A JSON5 string representing the value.

### Node.js `require()` JSON5 files
When using Node.js, you can `require()` JSON5 files by adding the following
statement.

```js
require('json5/lib/register')
```

Then you can load a JSON5 file with a Node.js `require()` statement. For
example:

```js
const config = require('./config.json5')
```

## CLI
Since JSON is more widely used than JSON5, this package includes a CLI for
converting JSON5 to JSON and for validating the syntax of JSON5 documents.

### Installation
```sh
npm install --global json5
```

### Usage
```sh
json5 [options] <file>
```

If `<file>` is not provided, then STDIN is used.

#### Options:
- `-s`, `--space`: The number of spaces to indent or `t` for tabs
- `-o`, `--out-file [file]`: Output to the specified file, otherwise STDOUT
- `-v`, `--validate`: Validate JSON5 but do not output JSON
- `-V`, `--version`: Output the version number
- `-h`, `--help`: Output usage information

## Contributing
### Development
```sh
git clone https://github.com/json5/json5
cd json5
npm install
```

When contributing code, please write relevant tests and run `npm test` and `npm
run lint` before submitting pull requests. Please use an editor that supports
[EditorConfig](http://editorconfig.org/).

### Issues
To report bugs or request features regarding the JSON5 data format, please
submit an issue to the [official specification
repository](https://github.com/json5/json5-spec).

To report bugs or request features regarding the JavaScript implementation of
JSON5, please submit an issue to this repository.

## License
MIT. See [LICENSE.md](./LICENSE.md) for details.

## Credits
[Assem Kishore](https://github.com/aseemk) founded this project.

[Michael Bolin](http://bolinfest.com/) independently arrived at and published
some of these same ideas with awesome explanations and detail. Recommended
reading: [Suggested Improvements to JSON](http://bolinfest.com/essays/json.html)

[Douglas Crockford](http://www.crockford.com/) of course designed and built
JSON, but his state machine diagrams on the [JSON website](http://json.org/), as
cheesy as it may sound, gave us motivation and confidence that building a new
parser to implement these ideas was within reach! The original
implementation of JSON5 was also modeled directly off of Doug’s open-source
[json_parse.js] parser. We’re grateful for that clean and well-documented
code.

[json_parse.js]:
https://github.com/douglascrockford/JSON-js/blob/03157639c7a7cddd2e9f032537f346f1a87c0f6d/json_parse.js

[Max Nanasy](https://github.com/MaxNanasy) has been an early and prolific
supporter, contributing multiple patches and ideas.

[Andrew Eisenberg](https://github.com/aeisenberg) contributed the original
`stringify` method.

[Jordan Tucker](https://github.com/jordanbtucker) has aligned JSON5 more closely
with ES5, wrote the official JSON5 specification, completely rewrote the
codebase from the ground up, and is actively maintaining this project.
# SQLite compiled to javascript
[![Build Status](https://travis-ci.org/kripken/sql.js.svg?branch=master)](http://travis-ci.org/kripken/sql.js) [![CDNJS version](https://img.shields.io/cdnjs/v/sql.js.svg)](https://cdnjs.com/libraries/sql.js)

For the impatients, try the demo here: http://kripken.github.io/sql.js/examples/GUI

*sql.js* is a port of [SQLite](http://sqlite.org/about.html) to Webassembly, by compiling the SQLite C code with [Emscripten](http://kripken.github.io/emscripten-site/docs/introducing_emscripten/about_emscripten.html). It uses a [virtual database file stored in memory](https://kripken.github.io/emscripten-site/docs/porting/files/file_systems_overview.html), and thus **doesn't persist the changes** made to the database. However, it allows you to **import** any existing sqlite file, and to **export** the created database as a [javascript typed array](https://developer.mozilla.org/en-US/docs/Web/JavaScript/Typed_arrays).

There are no C bindings or node-gyp compilation here, sql.js is a simple javascript file, that can be used like any traditional javascript library. If you are building a native application in javascript (using Electron for instance), or are working in node.js, you will likely prefer to use [a native binding of SQLite to javascript](https://www.npmjs.com/package/sqlite3).

SQLite is public domain, sql.js is MIT licensed.

Sql.js predates WebAssembly, and thus started as an [asm.js](https://en.wikipedia.org/wiki/Asm.js) project. It still supports asm.js for backwards compatability.

## Version of binaries
Sql.js was last built with:
Emscripten version 1.38.30 (2019-04-16) [Release History](https://emscripten.org/docs/introducing_emscripten/release_notes.html)
SqlLite version: 3.28.0 (2019-04-16) [Release History](https://www.sqlite.org/changes.html)

## Documentation
A [full documentation](http://kripken.github.io/sql.js/documentation/#http://kripken.github.io/sql.js/documentation/class/Database.html) generated from comments inside the source code, is available.

## Usage

```javascript
var initSqlJs = require('sql-wasm.js');
// or if you are in a browser:
//var initSqlJs = window.initSqlJs;

initSqlJs().then(function(SQL){

  // Create a database
  var db = new SQL.Database();
  // NOTE: You can also use new SQL.Database(data) where
  // data is an Uint8Array representing an SQLite database file

  // Execute some sql
  sqlstr = "CREATE TABLE hello (a int, b char);";
  sqlstr += "INSERT INTO hello VALUES (0, 'hello');"
  sqlstr += "INSERT INTO hello VALUES (1, 'world');"
  db.run(sqlstr); // Run the query without returning anything

  var res = db.exec("SELECT * FROM hello");
  /*
  [
    {columns:['a','b'], values:[[0,'hello'],[1,'world']]}
  ]
  */

  // Prepare an sql statement
  var stmt = db.prepare("SELECT * FROM hello WHERE a=:aval AND b=:bval");

  // Bind values to the parameters and fetch the results of the query
  var result = stmt.getAsObject({':aval' : 1, ':bval' : 'world'});
  console.log(result); // Will print {a:1, b:'world'}

  // Bind other values
  stmt.bind([0, 'hello']);
  while (stmt.step()) console.log(stmt.get()); // Will print [0, 'hello']

  // You can also use javascript functions inside your SQL code
  // Create the js function you need
  function add(a, b) {return a+b;}
  // Specifies the SQL function's name, the number of it's arguments, and the js function to use
  db.create_function("add_js", add);
  // Run a query in which the function is used
  db.run("INSERT INTO hello VALUES (add_js(7, 3), add_js('Hello ', 'world'));"); // Inserts 10 and 'Hello world'

  // free the memory used by the statement
  stmt.free();
  // You can not use your statement anymore once it has been freed.
  // But not freeing your statements causes memory leaks. You don't want that.

  // Export the database to an Uint8Array containing the SQLite database file
  var binaryArray = db.export();
});

```

## Demo
There are a few examples [available here](https://kripken.github.io/sql.js/index.html). The most full-featured is the [Sqlite Interpreter](https://kripken.github.io/sql.js/examples/GUI/index.html).

## Examples
The test files provide up to date example of the use of the api.
### Inside the browser
#### Example **HTML** file:
```html
<meta charset="utf8" />
<html>
  <script src='/dist/sql-wasm.js'></script>
  <script>
    config = {
      locateFile: url => `/dist/${filename}` 
    }
    // The `initSqlJs` function is globally provided by all of the main dist files if loaded in the browser.
    // We must specify this locateFile function if we are loading a wasm file from anywhere other than the current html page's folder.
    initSqlJs(config).then(function(SQL){
      //Create the database
      var db = new SQL.Database();
      // Run a query without reading the results
      db.run("CREATE TABLE test (col1, col2);");
      // Insert two rows: (1,111) and (2,222)
      db.run("INSERT INTO test VALUES (?,?), (?,?)", [1,111,2,222]);
  
      // Prepare a statement
      var stmt = db.prepare("SELECT * FROM test WHERE col1 BETWEEN $start AND $end");
      stmt.getAsObject({$start:1, $end:1}); // {col1:1, col2:111}
  
      // Bind new values
      stmt.bind({$start:1, $end:2});
      while(stmt.step()) { //
        var row = stmt.getAsObject();
        console.log('Here is a row: ' + JSON.stringify(row));
      }
    });
  </script>
  <body>
    Output is in Javscript console
  </body>
</html>
```

#### Creating a database from a file choosen by the user
`SQL.Database` constructor takes an array of integer representing a database file as an optional parameter.
The following code uses an HTML input as the source for loading a database:
```javascript
dbFileElm.onchange = () => {
  var f = dbFileElm.files[0];
  var r = new FileReader();
  r.onload = function() {
    var Uints = new Uint8Array(r.result);
    db = new SQL.Database(Uints);
  }
  r.readAsArrayBuffer(f);
}
```
See : http://kripken.github.io/sql.js/examples/GUI/gui.js

#### Loading a database from a server

```javascript
var xhr = new XMLHttpRequest();
// For example: https://github.com/lerocha/chinook-database/raw/master/ChinookDatabase/DataSources/Chinook_Sqlite.sqlite
xhr.open('GET', '/path/to/database.sqlite', true);
xhr.responseType = 'arraybuffer';

xhr.onload = e => {
  var uInt8Array = new Uint8Array(this.response);
  var db = new SQL.Database(uInt8Array);
  var contents = db.exec("SELECT * FROM my_table");
  // contents is now [{columns:['col1','col2',...], values:[[first row], [second row], ...]}]
};
xhr.send();
```
See: https://github.com/kripken/sql.js/wiki/Load-a-database-from-the-server


### Use from node.js

`sql.js` is [hosted on npm](https://www.npmjs.org/package/sql.js). To install it, you can simply run `npm install sql.js`.
Alternatively, you can simply download `sql-wasm.js` and `sql-wasm.wasm`, from the download link below.

#### read a database from the disk:
```javascript
var fs = require('fs');
var initSqlJs = require('sql-wasm.js');
var filebuffer = fs.readFileSync('test.sqlite');
 
initSqlJs().then(function(SQL){
  // Load the db
  var db = new SQL.Database(filebuffer);
});

```

#### write a database to the disk
You need to convert the result of `db.export` to a buffer
```javascript
var fs = require("fs");
// [...] (create the database)
var data = db.export();
var buffer = new Buffer(data);
fs.writeFileSync("filename.sqlite", buffer);
```

See : https://github.com/kripken/sql.js/blob/master/test/test_node_file.js

### Use as web worker
If you don't want to run CPU-intensive SQL queries in your main application thread,
you can use the *more limited* WebWorker API.

You will need to download [dist/worker.sql-wasm.js](dist/worker.sql-wasm.js) [dist/worker.sql-wasm.wasm](dist/worker.sql-wasm.wasm).

Example:
```html
<script>
  var worker = new Worker("/dist/worker.sql-wasm.js");
  worker.onmessage = () => {
    console.log("Database opened");
    worker.onmessage = event => {
      console.log(event.data); // The result of the query
    };
	
    worker.postMessage({
      id: 2,
      action: 'exec',
      sql: 'SELECT * FROM test'
    });
  };

  worker.onerror = e => console.log("Worker error: ", e);
  worker.postMessage({
    id:1,
    action:'open',
    buffer:buf, /*Optional. An ArrayBuffer representing an SQLite Database file*/
  });
</script>
```

See [examples/GUI/gui.js](examples/GUI/gui.js) for a full working example.

## Flavors/versions Targets/Downloads

This library includes both WebAssembly and asm.js versions of Sqlite. (WebAssembly is the newer, preferred way to compile to Javascript, and has superceded asm.js. It produces smaller, faster code.) Asm.js versions are included for compatibility.

## Upgrading from 0.x to 1.x

Version 1.0 of sql.js must be loaded asynchronously, whereas asm.js was able to be loaded synchronously. 

So in the past, you would:
```html
<script src='js/sql.js'></script>
<script>
  var db = new SQL.Database();
  //...
</script>
```
or:
```javascript
var SQL = require('sql.js');
var db = new QL.Database();
//...
```

Version 1.x:
```html
<script src='dist/sql-wasm.js'></script>
<script>
  initSqlJs({ locateFile: filename => `/dist/${filename}` }).then(function(SQL){
    var db = new SQL.Database();
    //...
  });
</script>
```
or:
```javascript
var initSqlJs = require('sql-wasm.js');
initSqlJs().then(function(SQL){
  var db = new SQL.Database();
  //...
});
```

`NOTHING` is now a reserved word in SQLite, whereas previously it was not. This could cause errors like `Error: near "nothing": syntax error`

### Downloading/Using: ###
Although asm.js files were distributed as a single Javascript file, WebAssembly libraries are most efficiently distributed as a pair of files, the `.js`  loader and the `.wasm` file, like [dist/sql-wasm.js]([dist/sql-wasm.js]) and [dist/sql-wasm.wasm]([dist/sql-wasm.wasm]). The `.js` file is reponsible for wrapping/loading the `.wasm` file. 




## Versions of sql.js included in `dist/`
 - `sql-wasm.js` : The Web Assembly version of Sql.js. Minified and suitable for production. Use this. If you use this, you will need to include/ship `sql-wasm.wasm` as well.
 - `sql-wasm-debug.js` : The Web Assembly, Debug version of Sql.js. Larger, with assertions turned on. Useful for local development. You will need to include/ship `sql-wasm-debug.wasm` if you use this.
 - `sql-asm.js` : The older asm.js version of Sql.js. Slower and larger. Provided for compatiblity reasons.
 - `sql-asm-memory-growth.js` : Asm.js doesn't allow for memory to grow by default, because it is slower and de-optimizes. If you are using sql-asm.js and you see this error (`Cannot enlarge memory arrays`), use this file.
 - `sql-asm-debug.js` : The _Debug_ asm.js version of Sql.js. Use this for local development.
 - `worker.*` - Web Worker versions of the above libraries. More limited API. See [examples/GUI/gui.js](examples/GUI/gui.js) for a good example of this.

## Compiling

- Install the EMSDK, [as described here](https://kripken.github.io/emscripten-site/docs/getting_started/downloads.html)
- Run `npm run rebuild`

# GeoPackage JS

GeoPackage JS is an implementation of the OGC GeoPackage spec.  This library works in both the browser and Node 4+.

### Demo ###
[GeoPackage JS Demo Page](http://ngageoint.github.io/geopackage-js/)

Cloning this repository and opening the docs/index.html in your browser will run the demo locally.

### Installation ###

[![Build Status](https://travis-ci.org/ngageoint/geopackage-js.svg?branch=master)](https://travis-ci.org/ngageoint/geopackage-js)
[![NPM](https://img.shields.io/npm/v/@ngageoint/geopackage.svg)](https://www.npmjs.com/package/@ngageoint/geopackage)
[![Coverage Status](https://coveralls.io/repos/github/ngageoint/geopackage-js/badge.svg)](https://coveralls.io/github/ngageoint/geopackage-js)

```sh
$ npm install @ngageoint/geopackage
```

#### GeoPackage JS Library ####

The [GeoPackage Libraries](http://ngageoint.github.io/GeoPackage/) were developed at the [National Geospatial-Intelligence Agency (NGA)](http://www.nga.mil/) in collaboration with [BIT Systems](http://www.bit-sys.com/). The government has "unlimited rights" and is releasing this software to increase the impact of government investments by providing developers with the opportunity to take things in new directions. The software use, modification, and distribution rights are stipulated within the [MIT license](http://choosealicense.com/licenses/mit/).

### Pull Requests ###
If you'd like to contribute to this project, please make a pull request. We'll review the pull request and discuss the changes. All pull request contributions to this project will be released under the MIT license.

Software source code previously released under an open source license and then modified by NGA staff is considered a "joint work" (see 17 USC § 101); it is partially copyrighted, partially public domain, and as a whole is protected by the copyrights of the non-government authors and must be released according to the terms of the original open source license.

### About ###

[GeoPackage JS](https://github.com/ngageoint/geopackage-js) is a [GeoPackage Library](http://ngageoint.github.io/GeoPackage/) JavaScript implementation of the Open Geospatial Consortium [GeoPackage](http://www.geopackage.org/) [spec](http://www.geopackage.org/spec/).  It is listed as an [OGC GeoPackage Implementation](http://www.geopackage.org/#implementations_nga) by the National Geospatial-Intelligence Agency.

The GeoPackage JavaScript library currently provides the ability to read GeoPackage files.  This library works both in the browser and in Node.  In the browser tiles are rendered using HTML5 Canvas and GeoPackages are read using [sql.js](https://github.com/kripken/sql.js/).  In Node tiles are rendered  [PureImage](https://github.com/joshmarinacci/node-pureimage) and GeoPackages are read using [node-sqlite3](https://github.com/mapbox/node-sqlite3).

### Changelog

##### 2.1.0

- Implementation of the Feature Style Extension and Contents ID Extension

##### 2.0.8

- Checks for Electron when returning a tile creator

##### 2.0

- All new API utilizing Promises

##### 1.1.4

- Adds a method to retrieve tiles in EPSG:4326

##### 1.1.3

- Fixes issue #115

##### 1.1.2

- fix case where GeoPackage Zoom does not correspond to the web map zoom

##### 1.1.1

- fix more instances of proj4 bug for react
- fixed tile generation for images with different x and y pixel densities

##### 1.1.0

- accept pull request adding support for react
- fix bug with projected tiles that spanned the date line

##### 1.0.25

- ensure we use proj4 2.4.3 instead of 2.4.4

##### 1.0.22

- Fixed bug where querying for indexed features only returned the geometry instead of the entire feature

##### 1.0.19

- Remove dependency on Lwip

### Usage ###

View examples using [Bower](https://github.com/ngageoint/geopackage-js/tree/master/docs/bower) and [Browserify](https://github.com/ngageoint/geopackage-js/tree/master/docs)

View the latest [docs](http://ngageoint.github.io/geopackage-js/jsdoc/module-geoPackage-GeoPackage.html) (currently being updated).

#### Browser Usage ####
```javascript

// attach this method to a file input onchange event
window.loadGeoPackage = function(files) {
  var f = files[0];
  var r = new FileReader();
  r.onload = function() {
    var array = new Uint8Array(r.result);
    loadByteArray(array);
  }
  r.readAsArrayBuffer(f);
}

function loadByteArray(array, callback) {
  var db = new SQL.Database(array);
  GeoPackageConnection.connectWithDatabase(db, function(err, connection) {
    var geoPackage = new GeoPackage('', '', connection);

    // Now you can operate on the GeoPackage

    // get the tile table names
    geoPackage.getTileTables(function(err, tileTableNames) {
      // tileTableNames is an array of all tile table names

      // get the info for the first table
      geoPackage.getTileDaoWithTableName(tileTableNames[0], function(err, tileDao) {
        geoPackage.getInfoForTable(tileDao, function(err, info) {
          // do something with the tile table info
        });

        // draw a tile into a canvas for an XYZ tile
        var canvas = canvasFromSomewhere;
        var gpr = new GeoPackageTileRetriever(tileDao, 256, 256);
        var x = 0;
        var y = 0;
        var zoom = 0;

        console.time('Draw tile ' + x + ', ' + y + ' zoom: ' + zoom);
        gpr.drawTileIn(x, y, zoom, canvas, function() {
          console.timeEnd('Draw tile ' + x + ', ' + y + ' zoom: ' + zoom);
        });

        // or get a tile base64 data URL for an XYZ tile
        gpr.getTile(x, y, zoom, function(err, tileBase64DataURL) {
          console.log('got the base64 data url');
        });

        // or get a tile from a GeoPackage tile column and tile row
        tileDao.queryForTile(tileColumn, tileRow, zoom, function(err, tile) {
          var tileData = tile.getTileData();  // the raw bytes from the GeoPackage
        });

      });
    });

    // get the feature table names
    geoPackage.getFeatureTables(function(err, featureTableNames) {
      // featureTableNames is an array of all feature table names

      // get the info for the first table
      geoPackage.getFeatureDaoWithTableName(featureTableNames[0], function(err, featureDao) {
        geoPackage.getInfoForTable(featureDao, function(err, info) {
          // do something with the feature table info
        });

        // query for all features
        featureDao.queryForEach(function(err, row, rowDone) {
          var feature = featureDao.getFeatureRow(row);
          var geometry = currentRow.getGeometry();
          if (geometry) {
            var geom = geometry.geometry;
            var geoJson = geometry.geometry.toGeoJSON();

            geoJson.properties = {};
            for (var key in feature.values) {
              if(feature.values.hasOwnProperty(key) && key != feature.getGeometryColumn().name) {
                var column = info.columnMap[key];
                geoJson.properties[column.displayName] = currentRow.values[key];
              }
            }
          }
          rowDone();
        });
      });
    });
  });
}

```

#### NodeJS Usage ####

```javascript
var GeoPackageAPI = require('@ngageoint/geopackage')
  , GeoPackageManager = GeoPackageAPI.GeoPackageManager
  , GeoPackageConnection = GeoPackageAPI.GeoPackageConnection
  , GeoPackageTileRetriever = GeoPackageAPI.GeoPackageTileRetriever;

GeoPackageAPI.open(filename, function(err, geoPackage) {

  // Now you can operate on the GeoPackage

  // get the tile table names
  geoPackage.getTileTables(function(err, tileTableNames) {
    // tileTableNames is an array of all tile table names

    // get the info for the first table
    geoPackage.getTileDaoWithTableName(tileTableNames[0], function(err, tileDao) {
      geoPackage.getInfoForTable(tileDao, function(err, info) {
        // do something with the tile table info
      });

      // draw a tile into a canvas for an XYZ tile
      var canvas = canvasFromSomewhere;
      var gpr = new GeoPackageTileRetriever(tileDao, 256, 256);
      var x = 0;
      var y = 0;
      var zoom = 0;

      console.time('Draw tile ' + x + ', ' + y + ' zoom: ' + zoom);
      gpr.drawTileIn(x, y, zoom, canvas, function() {
        console.timeEnd('Draw tile ' + x + ', ' + y + ' zoom: ' + zoom);
      });

      // or get a tile base64 data URL for an XYZ tile
      gpr.getTile(x, y, zoom, function(err, tileBase64DataURL) {
        console.log('got the base64 data url');
      });

      // or get a tile from a GeoPackage tile column and tile row
      tileDao.queryForTile(tileColumn, tileRow, zoom, function(err, tile) {
        var tileData = tile.getTileData();  // the raw bytes from the GeoPackage
      });

    });
  });

  // get the feature table names
  geoPackage.getFeatureTables(function(err, featureTableNames) {
    // featureTableNames is an array of all feature table names

    // get the info for the first table
    geoPackage.getFeatureDaoWithTableName(featureTableNames[0], function(err, featureDao) {
      geoPackage.getInfoForTable(featureDao, function(err, info) {
        // do something with the feature table info
      });

      // query for all features
      featureDao.queryForEach(function(err, row, rowDone) {
        var feature = featureDao.getFeatureRow(row);
        var geometry = currentRow.getGeometry();
        if (geometry) {
          var geom = geometry.geometry;
          var geoJson = geometry.geometry.toGeoJSON();

          geoJson.properties = {};
          for (var key in feature.values) {
            if(feature.values.hasOwnProperty(key) && key != feature.getGeometryColumn().name) {
              var column = info.columnMap[key];
              geoJson.properties[column.displayName] = currentRow.values[key];
            }
          }
        }
        rowDone();
      });
    });
  });
});

```
# PROJ4JS [![Build Status](https://travis-ci.org/proj4js/proj4js.svg)](https://travis-ci.org/proj4js/proj4js)

Proj4js is a JavaScript library to transform point coordinates from one coordinate system to another, including datum transformations.
Originally a port of [PROJ](https://proj.org/) ([then known as PROJ.4](https://proj.org/faq.html#what-happened-to-proj-4)) and GCTCP C ([Archive](https://web.archive.org/web/20130523091752/http://edcftp.cr.usgs.gov/pub/software/gctpc/)) it is
a part of the [MetaCRS](https://trac.osgeo.org/metacrs/wiki) group of projects.

## Installing

Depending on your preferences

```bash
npm install proj4
bower install proj4
component install proj4js/proj4js
```

or just manually grab the file `proj4.js` from the [latest release](https://github.com/proj4js/proj4js/releases)'s `dist/` folder.

If you do not want to download anything, Proj4js is also hosted on [cdnjs](https://www.cdnjs.com/libraries/proj4js) for direct use in your browser applications.

## Using

The basic signature is:

```javascript
proj4(fromProjection[, toProjection, coordinates])
```

Projections can be proj or wkt strings.

Coordinates may an object of the form `{x:x,y:y}` or an array of the form `[x,y]`.

When all 3 arguments  are given, the result is that the coordinates are transformed from projection1 to projection 2. And returned in the same format that they were given in.

```javascript
var firstProjection = 'PROJCS["NAD83 / Massachusetts Mainland",GEOGCS["NAD83",DATUM["North_American_Datum_1983",SPHEROID["GRS 1980",6378137,298.257222101,AUTHORITY["EPSG","7019"]],AUTHORITY["EPSG","6269"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4269"]],UNIT["metre",1,AUTHORITY["EPSG","9001"]],PROJECTION["Lambert_Conformal_Conic_2SP"],PARAMETER["standard_parallel_1",42.68333333333333],PARAMETER["standard_parallel_2",41.71666666666667],PARAMETER["latitude_of_origin",41],PARAMETER["central_meridian",-71.5],PARAMETER["false_easting",200000],PARAMETER["false_northing",750000],AUTHORITY["EPSG","26986"],AXIS["X",EAST],AXIS["Y",NORTH]]';
var secondProjection = "+proj=gnom +lat_0=90 +lon_0=0 +x_0=6300000 +y_0=6300000 +ellps=WGS84 +datum=WGS84 +units=m +no_defs";
//I'm not going to redefine those two in latter examples.
proj4(firstProjection,secondProjection,[2,5]);
// [-2690666.2977344505, 3662659.885459918]
```

If only 1 projection is given then it is assumed that it is being projected *from* WGS84 (fromProjection is WGS84).

```javascript
proj4(firstProjection,[-71,41]);
// [242075.00535055372, 750123.32090043]
```

If no coordinates are given an object with two methods is returned, its methods are `forward` which projects from the first projection to the second and `inverse` which projects from the second to the first.

```javascript
proj4(firstProjection,secondProjection).forward([2,5]);
// [-2690666.2977344505, 3662659.885459918]
proj4(secondProjection,firstProjection).inverse([2,5]);
// [-2690666.2977344505, 3662659.885459918]
```

And as above if only one projection is given, it's assumed to be coming from wgs84:

```javascript
proj4(firstProjection).forward([-71,41]);
// [242075.00535055372, 750123.32090043]
proj4(firstProjection).inverse([242075.00535055372, 750123.32090043]);
//[-71, 40.99999999999986]
//the floating points to answer your question
```

## Named Projections

If you prefer to define a projection as a string and reference it that way, you may use the proj4.defs method which can be called 2 ways, with a name and projection:

```js
proj4.defs('WGS84', "+title=WGS 84 (long/lat) +proj=longlat +ellps=WGS84 +datum=WGS84 +units=degrees");
```

or with an array

```js
proj4.defs([
  [
    'EPSG:4326',
    '+title=WGS 84 (long/lat) +proj=longlat +ellps=WGS84 +datum=WGS84 +units=degrees'],
  [
    'EPSG:4269',
    '+title=NAD83 (long/lat) +proj=longlat +a=6378137.0 +b=6356752.31414036 +ellps=GRS80 +datum=NAD83 +units=degrees'
  ]
]);
```

you can then do

```js
proj4('EPSG:4326');
```

instead of writing out the whole proj definition, by default proj4 has the following projections predefined:

- 'EPSG:4326', which has the following alias
    - 'WGS84'
- 'EPSG:4269'
- 'EPSG:3857', which has the following aliases
    - 'EPSG:3785'
    - 'GOOGLE'
    - 'EPSG:900913'
    - 'EPSG:102113'

Defined projections can also be accessed through the proj4.defs function (`proj4.defs('EPSG:4326')`).

proj4.defs can also be used to define a named alias:

```javascript
proj4.defs('urn:x-ogc:def:crs:EPSG:4326', proj4.defs('EPSG:4326'));
```

## TypeScript

TypeScript implementation was added to the [DefinitelyTyped repository](https://github.com/DefinitelyTyped/DefinitelyTyped).

```bash
$ npm install --save @types/proj4
```

## Developing
To set up build tools make sure you have node and grunt-cli installed and then run `npm install`.

To do the complete build and browser tests run:

```bash
node_modules/.bin/grunt
```

To run node tests run:

```bash
npm test
```

To run node tests with coverage run:

```bash
npm test --coverage
```

To create a build with only default projections (latlon and Mercator) run:

```bash
node_modules/.bin/grunt build
```

To create a build with only custom projections include a comma separated list of projections codes (the file name in 'lib/projections' without the '.js') after a colon, e.g.:

```bash
node_modules/.bin/grunt build:tmerc
#includes transverse Mercator
node_modules/.bin/grunt build:lcc
#includes lambert conformal conic
node_modules/.bin/grunt build:omerc,moll
#includes oblique Mercator and Mollweide
```
