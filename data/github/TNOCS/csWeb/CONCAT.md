[![Codacy Badge](https://api.codacy.com/project/badge/Grade/c1bf4babb31d4e7196dbda21088647d2)](https://www.codacy.com/app/erikvullings/csWeb?utm_source=github.com&utm_medium=referral&utm_content=TNOCS/csWeb&utm_campaign=badger)
[![Build Status](https://travis-ci.org/TNOCS/csWeb.svg?branch=master)](https://travis-ci.org/TNOCS/csWeb)
[![Stories in Ready](https://badge.waffle.io/tnocs/csweb.png?label=ready&title=Ready)](https://waffle.io/tnocs/csweb)
[![bitHound Score](https://www.bithound.io/github/TNOCS/csWeb/badges/score.svg)](https://www.bithound.io/github/TNOCS/csWeb/layer-sources-renders)
[![Coverage Status](https://coveralls.io/repos/TNOCS/csWeb/badge.svg?branch=development)](https://coveralls.io/r/TNOCS/csWeb?branch=development)

# README #

**csWeb**, or the **Common Sense Web application**, is an intuitive open source web-based GIS application, providing casual users as well as business analysists and information manageners with a powerful tool to perform spatial analysis. It has a strong focus on usability and connectivity, be it connecting and sharing information with other users or connecting to services or calculation simulations and models. [LIVE DEMO](http://tnocs.github.io/csWeb/)

## Features
* Basic map interactions (zooming, geo-locating, selecting different base layers)
* Displaying geojson files
* Specifying how properties must be displayed (formatting, title, tooltips, etc.)
* Filtering on one or more properties
* Coloring icons and regions based on one or more properties
* Displaying data in a table, and allowing the users to download it
* Searching for a feature

## Technical overview

Technically, we use the following frameworks:
* Typescript for coding the application
* Angularjs as the MVC framework
* Bootstrap 3 (and fontawesome) for the CSS design style
* Leaflet for the 2D map layer, and Cesium for the 3D maps
* d3, dc, crossfilter for filtering and styling

### How do I get set up? ###

The application is written in Typescript, which compiles to regular JavaScript, and further uses Angularjs as the framework, Leaflet and Cesium for rendering maps, d3 and others. A detailed guide is provided [here](https://github.com/TNOCS/csWeb/wiki/Getting-started), a concise [installation checklist here](https://github.com/TNOCS/csWeb/wiki/Installation-checklist).

This repository consists of several project folders. The  most important ones are csComp, a library containing client side functionality, and csServerComp, for server side components. Both libraries generate a JavaScript file that can be used by the actual map application which you can find in the `example` folder.

#### Deployment instructions ####

Just copy the example folder to a public folder and open the public\index.html file in that folder.

### Using Docker containers
There are two types of docker containers built for CommonSense:
* [`tnocs/csWeb-demo`](https://hub.docker.com/r/tnocs/csweb-demo/)
* [`tnocs/csWeb-dev`](https://hub.docker.com/r/tnocs/csweb-dev/)

#### `tnocs/csWeb-demo`
This container runs whole csWeb application. To run it and access csWeb in the web browser at `<port>` run:
```sh
# replace <port> with the port number you want to access csWeb at
docker run -d -p 3002:<port> tnocs/csWeb-demo
```
If your're using docker running in your system (Linux), application should be avaiilable at `localhost`,
otherwise when using `docker-machine`, you can ask about ip address with `docker-machine ip default` (in case your vm with docker is named default).

#### `tnocs/csWeb-dev`
This container is meant to run csWeb build in the docker container istread of local machine. csWeb local repository will be mounted inside the container and the build process will happen inside the container. This is usefull in at least two cases.
* to check if you understand all dependencies of your build process since only dependencies specified in the [Dockerfile](https://github.com/TNOCS/csWeb/blob/development-docker/docker-dev/Dockerfile) will be installed
* to avoid installing dependencies on local machin

Run container from within csWeb directory:
```sh
docker run -it --rm -p 3002:<port> -v $PWD:/app/ tno/csWeb-dev /bin/bash
```

### Contribution guidelines ###

* Writing tests
* Code review
* Other guidelines

#### Testing ####



```
karma start test/karma.conf.js
```

### Who do I talk to? ###

* Repo owner or admin
﻿# csServerComp


# Type of rules

## Conditions

### Time based

1. After a certain time
2. At a certain time

### Event based

An event is either a key/namespace of the event, optionally a subject, and a value

1. When an event is raised, optionally with a certain subject.

### Property based

1. When a certain property

## Actions

### Sending a message (feature update)

### Activating or deactivating rules


## Examples

* Time: Time progresses and messages / feature updates are sent.

"_rules": [{
    "id": "rule1",
    "desc": "Add the item after a certain amount of time (after activation)",
    "isActive": true,
    "actions": [{ "add", 20 }] // add the feature after 60s (that the condition becomes true)
}]

* Request: A question is asked about the social security number. The system responds, after a delay, with an answer.

"_rules": [{
    "actions": [{ "add", 60 }] // add the feature after 60s (that the condition becomes true)
}, {
    "conditions": [{ "property", "bsn" }, { "property", "isAnswered", false }, { "property", "assigned" }]
    "actions": [{ "set", "bsn", 123, 30 }, { "set", "isAnswered", true, 30 }] // set property to value after delay
}]

event: { "name": "bsn", "subject": "question", "isAnswered": false }
During the evaluation of the event, the feature.properties["isAnswered"] is compared to false.

action (context = feature): activateRule(feature, "isAnswered": true, "bsn": "123", delay = 60s);

* Request: The walking route is requested. In response, a route (polyline) is provided.

"_rules": [{
    "actions": [{ "add", 20 }] // add the feature after 60s (that the condition becomes true)
}, {
    "conditions": [{ "property", "bsn" }, { "property", "isAnswered", false }, { "property", "assigned" }]
    "actions": [{ "set", "bsn", 123, 30 }, { "set", "isAnswered", true, 30 }] // set property to value after delay
}]

conditions: [{ "property", "route" }, { "property", "isAnswered", false }, { "property", "assigned" }]
actions: [{ "add", 60 }] // add the feature after 60s (that the condition becomes true)

Do we evaluate the condition again before executing the action, e.g. when it has been recalled?.3..

.0.
0..

0..
0.20

* Request: An image is requested. In response, an image is provided.

* Action: A site or shop should be inspected by an actor. In response, the site/shop is visited by an actor and a result is sent (kid (not) found).

* Action: The wife / neighbor / school / wijkagent should be called. The result of this call is returned.
# csComp - Client side components

csComp contains the client side Angular components.

﻿The language switch directive allows you to switch the GUI language interactively. Assuming you have created
a translation. you can add it easilyby taking the following steps:

1. Go to [FAMFAMFAM](http://www.famfamfam.com/lab/icons/flags/) and download the flag you need.
2. Go to [Base64 encoder](http://www.base64-image.de/) to encode your flag's png image.
3. Add the XHTML Image string to the languagesProvider translations in app.ts.

    angular.module('csWebApp', [
            'pascalprecht.translate',
            'csWeb.languageSwitch'
        ])
        .config($translateProvider => {
            // TODO ADD YOUR LOCAL TRANSLATIONS HERE, OR ALTERNATIVELY, CHECK OUT 
            // http://angular-translate.github.io/docs/#/guide/12_asynchronous-loading
            // Translations.English.locale['MAP_LABEL'] = 'MY AWESOME MAP';
            $translateProvider.translations('en', Translations.English.locale);
            $translateProvider.translations('nl', Translations.Dutch.locale);
            $translateProvider.preferredLanguage('en');
        })
        .config(languagesProvider => {
            // Defines the GUI languages that you wish to use in your project.
            // They will be available through a popup menu.
            var languages = [];
            languages.push({ key: 'en', name: 'English'   , img: 'data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAABAAAAALCAIAAAD5gJpuAAAAGXRFWHRTb2Z0d2FyZQBBZG9iZSBJbWFnZVJlYWR5ccllPAAAAflJREFUeNpinDRzn5qN3uFDt16+YWBg+Pv339+KGN0rbVP+//2rW5tf0Hfy/2+mr99+yKpyOl3Ydt8njEWIn8f9zj639NC7j78eP//8739GVUUhNUNuhl8//ysKeZrJ/v7z10Zb2PTQTIY1XZO2Xmfad+f7XgkXxuUrVB6cjPVXef78JyMjA8PFuwyX7gAZj97+T2e9o3d4BWNp84K1NzubTjAB3fH0+fv6N3qP/ir9bW6ozNQCijB8/8zw/TuQ7r4/ndvN5mZgkpPXiis3Pv34+ZPh5t23//79Rwehof/9/NDEgMrOXHvJcrllgpoRN8PFOwy/fzP8+gUlgZI/f/5xcPj/69e/37//AUX+/mXRkN555gsOG2xt/5hZQMwF4r9///75++f3nz8nr75gSms82jfvQnT6zqvXPjC8e/srJQHo9P9fvwNtAHmG4f8zZ6dDc3bIyM2LTNlsbtfM9OPHH3FhtqUz3eXX9H+cOy9ZMB2o6t/Pn0DHMPz/b+2wXGTvPlPGFxdcD+mZyjP8+8MUE6sa7a/xo6Pykn1s4zdzIZ6///8zMGpKM2pKAB0jqy4UE7/msKat6Jw5mafrsxNtWZ6/fjvNLW29qv25pQd///n+5+/fxDDVbcc//P/zx/36m5Ub9zL8+7t66yEROcHK7q5bldMBAgwADcRBCuVLfoEAAAAASUVORK5CYII=' });
            languages.push({ key: 'nl', name: 'Nederlands', img: 'data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAABAAAAALCAIAAAD5gJpuAAAABGdBTUEAAK/INwWK6QAAABl0RVh0U29mdHdhcmUAQWRvYmUgSW1hZ2VSZWFkeXHJZTwAAAFXSURBVHjaYvzPgAD/UNlYEUAAkuTgCAAIBgJggq5VoAs1qM0vdzmMz362vezjokxPGimkEQ5WoAQEKuK71zwCCKyB4c//J8+BShn+/vv/+w/D399AEox+//8FJH/9/wUU+cUoKw20ASCAWBhEDf/LyDOw84BU//kDtgGI/oARmAHRDJQSFwVqAAggxo8fP/Ly8oKc9P8/AxjiAoyMjA8ePAAIIJZ///5BVIM0MOBWDpRlZPzz5w9AALH8gyvCbz7QBrCJAAHEyKDYX15r/+j1199//v35++/Xn7+///77DST/wMl/f4Dk378K4jx7O2cABBALw7NP77/+ev3xB0gOpOHfr99AdX9/gTVASKCGP//+8XCyMjC8AwggFoZfIHWSwpwQk4CW/AYjsKlA8u+ff////v33998/YPgBnQQQQIzAaGNg+AVGf5AYf5BE/oCjGEIyAQQYAGvKZ4C6+xXRAAAAAElFTkSuQmCC' });
            languages.push({ key: 'de', name: 'Deutsch'   , img: 'data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAABAAAAALCAIAAAD5gJpuAAAABGdBTUEAAK/INwWK6QAAABl0RVh0U29mdHdhcmUAQWRvYmUgSW1hZ2VSZWFkeXHJZTwAAAGzSURBVHjaYvTxcWb4+53h3z8GZpZff/79+v3n/7/fDAz/GHAAgABi+f37e3FxOZD1Dwz+/v3z9y+E/AMFv3//+Qumfv9et241QACxMDExAVWfOHkJJAEW/gUEP0EQDn78+AHE/gFOQJUAAcQiy8Ag8O+fLFj1n1+/QDp+/gQioK7fP378+vkDqOH39x9A/RJ/gE5lAAhAYhzcAACCQBDkgRXRjP034R0IaDTZTFZn0DItot37S94KLOINerEcI7aKHAHE8v/3r/9//zIA1f36/R+o4tevf1ANYNVA9P07RD9IJQMDQACxADHD3z8Ig4GMHz+AqqHagKp//fwLVA0U//v7LwMDQACx/LZiYFD7/5/53/+///79BqK/EMZ/UPACSYa/v/8DyX9A0oTxx2EGgABi+a/H8F/m339BoCoQ+g8kgRaCQvgPJJiBYmAuw39hxn+uDAABxMLwi+E/0PusRkwMvxhBGoDkH4b/v/+D2EDyz///QB1/QLb8+sP0lQEggFh+vGXYM2/SP6A2Zoaf30Ex/J+PgekHwz9gQDAz/P0FYrAyMfz7wcDAzPDtFwNAgAEAd3SIyRitX1gAAAAASUVORK5CYII=' });
            languagesProvider.setLanguages(languages);
        })

I've added the directive to the index page, but you can take your own pick.

    <!-- Language switch -->
    <language-switch class="navbar-form navbar-right" style="margin-top: 0; margin-bottom: 0"></language-switch>
    <!-- Search form --
= Grid-based heatmap algorithm =

== Introduction ==

Normally, a heatmap represents a map where the density of an item is depicted. So where we have many items of type x, the area is colored red (or orange, yellow green when the density decreases). However, in this case we wish to create a heatmap that represents the ideality of a location: based on certain criteria, color an area, thereby assuming that the presence of a certain feature on the map has a positive of negative influence on that area. For an example of a server-side solution, see geotrellis.io. 

The difference with a multi-criteria analysis (MCA) is that in an MCA, we compute the ideality of the item itself, and in an ideality map, we compute the ideality of an area on the map.   

== Algorithm ==

Compute a heatmap a uniform grid for the heatmap (heigth = width in meters of a cell, but not uniform in degrees lat/lon) on the map

For example, at 50 degrees latitude, one degree is:
lat 111229.03 meter (North-South)
lon  71695.73 meter (East-West)

At 57 degrees latitude, one degree is:
lat 111359.83 meter (North-South)
lon  60772.16 meter (East-West)

As you can see, there especially is a big difference in the width (East-West) of a grid cell when moving North. 

1. Determine the map extent
2. Determine the height in meter -> choose a grid cell size which covers the extent (max 100? rows)
3. For every selected HeatmapItem, compute a matrix [i j weighted_intensity], where i and j are indexes in a grid cell relative to the centre of the spot, and the intensity is:

	* for i=j=0 is the average intensity between the centre and the border of the grid 
	* for all other i,j is the intensity based on the distance between the centre and the centre of the current cell
	* the weighted intensity is the computed intensity * weight / scale
	* the scale is a relative factor: when scale is 5, it means that we need 5 spots at the ideal distance in a grid cell (we don't want to reach a maximum intensity when one item is at the ideal distance.)

4. Expand the extent by x meters in each direction (so we will include features just outside the border too), where x is determined by the HeatmapItem with the largest radius of influence.
5. For every feature within the expanded extent, add the intensity to the grid matrix (also [i j total_intensity])
6. Compute a color scale, mapping total_intensity (in the range [-1 1]) to a color [red white blue], where white reresents 0. In case the total_intensity exceeds the range, cap it to [-1 1]. Alternatively, first compute the max and min, and map them to blue and red respectively, where a value of 0 indicates transparent.  
6. Create a GeoJSON grid, where each height and width are equal to the grid cell size (this means that the longitude angle will be wider when moving North), and whose color is determined by the grid matrix's total_intensity value. In case the total_intensity is 0 (or a very small value), ignore the cell and don't add it to the GeoJSON grid.

When panning or zooming in/out, recalculate the grid, and replace the current grid when you are done.# angular-utils-pagination

This repo contains only the release version of the dirPagination module from my
[angularUtils](https://github.com/michaelbromley/angularUtils/tree/master/src/directives/pagination) repo. The sole purpose of this repo is to allow for convenient 
dependency management via Bower, npm and other package managers.

## Documentation

All documentation is also located in the [angularUtils project](https://github.com/michaelbromley/angularUtils/tree/master/src/directives/pagination).

## Issues

Please submit any issues to the [angularUtils issue tracker](https://github.com/michaelbromley/angularUtils/issues), not this one.

## Contribution

If you wish to contribute, then please fork the [angularUtils repo](https://github.com/michaelbromley/angularUtils). This repo contains
the tests and is set up to run Karma for unit testing.

## License 

MIT
# packaged angular-mocks

This repo is for distribution on `npm` and `bower`. The source for this module is in the
[main AngularJS repo](https://github.com/angular/angular.js/tree/master/src/ngMock).
Please file issues and pull requests against that repo.

## Install

You can install this package either with `npm` or with `bower`.

### npm

```shell
npm install angular-mocks
```

You can `require` ngMock modules:

```js
var angular = require('angular');
angular.module('myMod', [
  require('angular-animate'),
  require('angular-mocks/ngMock')
  require('angular-mocks/ngAnimateMock')
]);
```

### bower

```shell
bower install angular-mocks
```

The mocks are then available at `bower_components/angular-mocks/angular-mocks.js`.

## Documentation

Documentation is available on the
[AngularJS docs site](https://docs.angularjs.org/guide/unit-testing).

## License

The MIT License

Copyright (c) 2010-2015 Google, Inc. http://angularjs.org

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
