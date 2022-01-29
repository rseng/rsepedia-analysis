<div align="center">
    <img src="https://raw.githubusercontent.com/JuliaPlots/Makie.jl/master/assets/makie_logo_canvas.svg" alt="Makie.jl">
</div>

From the japanese word [_Maki-e_](https://en.wikipedia.org/wiki/Maki-e), which is a technique to sprinkle lacquer with gold and silver powder.
Data is the gold and silver of our age, so let's spread it out beautifully on the screen!

[Check out the documentation here!](http://makie.juliaplots.org/stable/)

[![][docs-stable-img]][docs-stable-url] [![][docs-master-img]][docs-master-url]

[gitlab-img]: https://gitlab.com/JuliaGPU/Makie.jl/badges/master/pipeline.svg
[gitlab-url]: https://gitlab.com/JuliaGPU/Makie.jl/pipelines
[docs-stable-img]: https://img.shields.io/badge/docs-stable-lightgrey.svg
[docs-stable-url]: http://makie.juliaplots.org/stable/
[docs-master-img]: https://img.shields.io/badge/docs-master-blue.svg
[docs-master-url]: http://makie.juliaplots.org/dev/

# Citing Makie

If you use Makie for a scientific publication, please cite [our JOSS paper](https://joss.theoj.org/papers/10.21105/joss.03349) the following way:

> Danisch & Krumbiegel, (2021). Makie.jl: Flexible high-performance data visualization for Julia. Journal of Open Source Software, 6(65), 3349, https://doi.org/10.21105/joss.03349

BibTeX entry:

```bib
@article{DanischKrumbiegel2021,
  doi = {10.21105/joss.03349},
  url = {https://doi.org/10.21105/joss.03349},
  year = {2021},
  publisher = {The Open Journal},
  volume = {6},
  number = {65},
  pages = {3349},
  author = {Simon Danisch and Julius Krumbiegel},
  title = {Makie.jl: Flexible high-performance data visualization for Julia},
  journal = {Journal of Open Source Software}
}
```

# Installation

Please consider using the backends directly. As explained in the documentation, they re-export all of Makie's functionality.
So, instead of installing Makie, just install e.g. GLMakie directly:
```julia
julia>]
pkg> add GLMakie
pkg> test GLMakie
```


Interactive example by [AlexisRenchon](https://github.com/AlexisRenchon):

![out](https://user-images.githubusercontent.com/1010467/81500379-2e8cfa80-92d2-11ea-884a-7069d401e5d0.gif)

Example from [InteractiveChaos.jl](https://github.com/JuliaDynamics/InteractiveChaos.jl)

[![interactive chaos](https://user-images.githubusercontent.com/1010467/81500069-ea005f80-92cf-11ea-81db-2b7bcbfea297.gif)
](https://github.com/JuliaDynamics/InteractiveChaos.jl)


You can follow Makie on [twitter](https://twitter.com/MakiePlots) to get the latest, outstanding examples:
[![image](https://user-images.githubusercontent.com/1010467/81500210-e7523a00-92d0-11ea-9849-1240f165e0f8.png)](https://twitter.com/MakiePlots)


## Sponsors

<img src="https://github.com/JuliaPlots/Makie.jl/blob/master/assets/BMBF_gefoerdert_2017_en.jpg?raw=true" width="300"/>
Förderkennzeichen: 01IS10S27, 2020
# News

## v0.16

#### Big Changes

- add ECDF plot [#1310](https://github.com/JuliaPlots/Makie.jl/pull/1310)
- add Order Independent Transparency to GLMakie [#1418](https://github.com/JuliaPlots/Makie.jl/pull/1418), [#1506](https://github.com/JuliaPlots/Makie.jl/pull/1506). This type of transparency is now used with `transpareny = true`. The old transparency handling is available with `transparency = false`.
- fix blurry text in GLMakie and WGLMakie [#1494](https://github.com/JuliaPlots/Makie.jl/pull/1494)
- A new experimental Backend for ray tracing got introduced: [RPRMakie](https://makie.juliaplots.org/stable/documentation/backends/rprmakie/)
- **Breaking** Remove `Node` alias [#1307](https://github.com/JuliaPlots/Makie.jl/pull/1307), [#1393](https://github.com/JuliaPlots/Makie.jl/pull/1393). To upgrade, simply replace all occurrences of `Node` with `Observable`
- **Breaking** clean up Scene type [#1192](https://github.com/JuliaPlots/Makie.jl/pull/1192), [#1393](https://github.com/JuliaPlots/Makie.jl/pull/1393). Long story short, Scene() doesn't create any axes or limits anymore. All keywords like `raw`, `show_axis` have been removed. A scene now always works like when using the deprecated `raw=true`. All the high level functionality like showing an axis and adding a 3d camera has been moved to `LScene`. See the new `Scene` tutorial for more info: https://makie.juliaplots.org/dev/tutorials/scenes/
- **lights got moved to scene** [lighting docs](https://makie.juliaplots.org/stable/documentation/lighting) and [RPRMakie examples](https://makie.juliaplots.org/stable/documentation/backends/rprmakie/)


#### Small Changes

- Added the `Cycled` type, which can be used to select the i-th value from the current cycler for a specific attribute. [#1248](https://github.com/JuliaPlots/Makie.jl/pull/1248)
- The plot function `scatterlines` now uses `color` as `markercolor` if `markercolor` is `automatic`. Also, cycling of the `color` attribute is enabled. [#1463](https://github.com/JuliaPlots/Makie.jl/pull/1463)
- Added the function `resize_to_layout!`, which allows to resize a `Figure` so that it contains its top `GridLayout` without additional whitespace or clipping. [#1438](https://github.com/JuliaPlots/Makie.jl/pull/1438)
- Cleanup lighting in 3D contours and isosurfaces [#1434](https://github.com/JuliaPlots/Makie.jl/pull/1434)
- Adjust attributes of volumeslices to follow the normal structure [#1404](https://github.com/JuliaPlots/Makie.jl/pull/1404). This allows you to adjust attributes like `colormap` without going through nested attributes.
- Add depth to 3D contours and isosurfaces [#1395](https://github.com/JuliaPlots/Makie.jl/pull/1395), [#1393](https://github.com/JuliaPlots/Makie.jl/pull/1393). This allows them to intersect correctly with other 3D objects.
- Restrict 3D scene camera to one scene [#1394](https://github.com/JuliaPlots/Makie.jl/pull/1394), [#1393](https://github.com/JuliaPlots/Makie.jl/pull/1393). This fixes issues with multiple scenes fighting over events consumed by the camera. You can select a scene by cleaning on it.
- add depth shift attribute for GLMakie and WGLMakie [#1382](https://github.com/JuliaPlots/Makie.jl/pull/1382), [#1393](https://github.com/JuliaPlots/Makie.jl/pull/1393). This can used to adjust render order similar to `overdraw`
- simplify automatic width computation in barplots [#1223](https://github.com/JuliaPlots/Makie.jl/pull/1223), [#1393](https://github.com/JuliaPlots/Makie.jl/pull/1393). If no `width` attribute is passed, the default width is computed as the minimum difference between consecutive `x` positions. Gap between bars are given by the (multiplicative) `gap` attribute. The actual bar width equals `width * (1 - gap)`.
- add logical expressions for `ispressed` [#1222](https://github.com/JuliaPlots/Makie.jl/pull/1222), [#1393](https://github.com/JuliaPlots/Makie.jl/pull/1393). This moves a lot of control over hotkeys towards the user. With these changes one can now set an hotkey to trigger on any or no key, collections of keys and logical combinations of keys (i.e. "A is pressed and B is not pressed").
- fix issues with Menu render order [#1411](https://github.com/JuliaPlots/Makie.jl/pull/1411)
- add label_rotation to barplot [#1401](https://github.com/JuliaPlots/Makie.jl/pull/1401)
- fix issue where `pixelcam!` does not remove controls from other cameras [#1504](https://github.com/JuliaPlots/Makie.jl/pull/1504)
- add conversion for offsetarrays [#1260](https://github.com/JuliaPlots/Makie.jl/pull/1260)
- `qqplot` `qqline` options are now `:identity`, `:fit`, `:fitrobust` and `:none` (the default) [#1563](https://github.com/JuliaPlots/Makie.jl/pull/1563). Fixed numeric error due to double computation of quantiles when fitting `qqline`. Deprecated `plot(q::QQPair)` method as it does not have enough information for correct `qqline` fit.

#### All other changes
Are collected [in this PR](https://github.com/JuliaPlots/Makie.jl/pull/1521) and in the [release notes](https://github.com/JuliaPlots/Makie.jl/releases/tag/v0.16.0).

## v0.15.3
- The functions `labelslidergrid!` and `labelslider!` now set fixed widths for the value column with a heuristic. It is possible now to pass `Formatting.format` format strings as format specifiers in addition to the previous functions.
- fix 2D arrow rotations in `streamplot` [#1352](https://github.com/JuliaPlots/Makie.jl/pull/1352)

## v0.15.2
- Reenabled Julia 1.3 support.
- Use [MathTexEngine v0.2](https://github.com/Kolaru/MathTeXEngine.jl/releases/tag/v0.2.0).
- Depend on new GeometryBasics, which changes all the Vec/Point/Quaternion/RGB/RGBA - f0 aliases to just f. For example, `Vec2f0` is changed to `Vec2f`. Old aliases are still exported, but deprecated and will be removed in the next breaking release. For more details and an upgrade script, visit [GeometryBasics#97](https://github.com/JuliaGeometry/GeometryBasics.jl/pull/97).
- Added `hspan!` and `vspan!` functions [#1264](https://github.com/JuliaPlots/Makie.jl/pull/1264).

## v0.15.1
- Switched documentation framework to Franklin.jl.
- Added a specialization for `volumeslices` to DataInspector.
- Fix [1 element `hist`](https://github.com/JuliaPlots/Makie.jl/pull/1238) and make it [easier to move `hist`](https://github.com/JuliaPlots/Makie.jl/pull/1150).

## v0.15.0

- `LaTeXString`s can now be used as input to `text` and therefore as labels for `Axis`, `Legend`, or other comparable objects. Mathematical expressions are typeset using [MathTeXEngine.jl](https://github.com/Kolaru/MathTeXEngine.jl) which offers a fast approximation of LaTeX typesetting. [#1022](https://github.com/JuliaPlots/Makie.jl/pull/1022)
- Added `Symlog10` and `pseudolog10` axis scales for log scale approximations that work with zero and negative values. [#1109](https://github.com/JuliaPlots/Makie.jl/pull/1109)
- Colorbar limits can now be passed as the attribute `colorrange` similar to plots. [#1066](https://github.com/JuliaPlots/Makie.jl/pull/1066)
- Added the option to pass three vectors to heatmaps and other plots using `SurfaceLike` conversion. [#1101](https://github.com/JuliaPlots/Makie.jl/pull/1101)
- Added `stairs` plot recipe. [#1086](https://github.com/JuliaPlots/Makie.jl/pull/1086)
- Removed `FigurePosition` and `FigureSubposition` types. Indexing into a `Figure` like `fig[1, 1]` now returns `GridPosition` and `GridSubposition` structs, which can be used in the same way as the types they replace. Because of an underlying change in `GridLayoutBase.jl`, it is now possible to do `Axis(gl[1, 1])` where `gl` is a `GridLayout` that is a sublayout of a `Figure`'s top layout. [#1075](https://github.com/JuliaPlots/Makie.jl/pull/1075)
- Bar plots and histograms have a new option for adding text labels. [#1069](https://github.com/JuliaPlots/Makie.jl/pull/1069)
- It is possible to specify one linewidth value per segment in `linesegments`. [#992](https://github.com/JuliaPlots/Makie.jl/pull/992)
- Added a new 3d camera that allows for better camera movements using keyboard and mouse. [#1024](https://github.com/JuliaPlots/Makie.jl/pull/1024)
- Fixed the application of scale transformations to `surface`. [#1070](https://github.com/JuliaPlots/Makie.jl/pull/1070)
- Added an option to set a custom callback function for the `RectangleZoom` axis interaction to enable other use cases than zooming. [#1104](https://github.com/JuliaPlots/Makie.jl/pull/1104)
- Fixed rendering of `heatmap`s with one or more reversed ranges in CairoMakie, as in `heatmap(1:10, 10:-1:1, rand(10, 10))`. [#1100](https://github.com/JuliaPlots/Makie.jl/pull/1100)
- fixed volume slice recipe and add docs for it [#1123](https://github.com/JuliaPlots/Makie.jl/pull/1123)
Makie Community Standards
=========================

The Makie community abides by the larger Julia Communities Code of Conduct (CoC). You can find that CoC listed here: https://julialang.org/community/standards/

If you have a conflict or concern that requires resolution, please contact the [Julia Community Stewards](https://julialang.org/community/stewards/).
# Contribution guidelines for Makie

## Issues

Before filing an issue please
- check that there are no similar existing issues already
- check that your versions are up to date
  - you can do this by forcing the latest version using `]add Makie@newest-version`, where `newest-version` can be found on the [releases page](https://github.com/JuliaPlots/Makie.jl/releases).

If you want to report a bug, include your version and system information, as well as stack traces with all relevant information.
If possible, condense your bug into the shortest example possible that the maintainers can replicate, a so called "minimal working example" or MWE.

If you want to report a visual regression, include an image showing the faulty behavior.

If you want to suggest a new feature, for example functionality that other plotting packages offer already, include supplementary material such as example images if possible, so it's clear what you are asking for.

## Pull requests

We always appreciate pull requests that fix bugs, add tests, improve documentation or add new features.
Please describe the intent of your pull request clearly and keep a clean commit history as far as possible.

For each feature you want to contribute, please file a separate PR to keep the complexity down and time to merge short.
Add PRs in draft mode if you want to discuss your approach first.

Please add tests for any new functionality that you want to add.
Makie uses both reference tests that check for visual regressions, and unit tests that check correctness of functions etc.
It is also appreciated if you add docstrings or documentation, and add an entry to the NEWS file.

### Tests

Please ensure locally that your feature works by running the tests.
To be able to run the tests, you have to `dev` the helper package `ReferenceTests` that is part of the Makie monorepo.
`ReferenceTests` is not a registered package, so you have to do `]dev path/to/ReferenceTests`.

After that you should be able to run the tests via the usual `]test` commands.

## Seeking Help

If you get stuck, here are some options to seek help:

- Use the REPL `?` help mode.
- Click this link to open a preformatted topic on the [Julia Discourse Page](https://discourse.julialang.org/new-topic?title=Makie%20-%20Your%20question%20here&category=domain/viz&tags=Makie&body=You%20can%20write%20your%20question%20in%20this%20space.%0A%0ABefore%20asking%2C%20please%20take%20a%20minute%20to%20make%20sure%20that%20you%20have%20installed%20the%20latest%20available%20versions%20and%20have%20looked%20at%20%5Bthe%20most%20recent%20documentation%5D(http%3A%2Fmakie.juliaplots.org%2Fstable%2F)%20%3Ainnocent%3A). If you do this manually, please use the category Domain/Visualization and tag questions with `Makie` to increase their visibility.
- For casual conversation about Makie and its development, have a look at the `#makie` channel in the [Julia Slack group](https://julialang.org/slack/). Please direct your usage questions to [Discourse](https://discourse.julialang.org/new-topic?title=Makie%20-%20Your%20question%20here&category=domain/viz&tags=Makie&body=You%20can%20write%20your%20question%20in%20this%20space.%0A%0ABefore%20asking%2C%20please%20take%20a%20minute%20to%20make%20sure%20that%20you%20have%20installed%20the%20latest%20available%20versions%20and%20have%20looked%20at%20%5Bthe%20most%20recent%20documentation%5D(http%3A%2Fmakie.juliaplots.org%2Fstable%2F)%20%3Ainnocent%3A) and not to Slack, to make questions and answers accessible to everybody.
# RPRMakie
# Description

Fixes # (issue)

Add a description of your PR here.

## Type of change

Delete options that do not apply:

- [ ] Bug fix (non-breaking change which fixes an issue)
- [ ] New feature (non-breaking change which adds functionality)
- [ ] Breaking change (fix or feature that would cause existing functionality to not work as expected)

## Checklist

- [ ] Added an entry in NEWS.md (for new features and breaking changes)
- [ ] Added or changed relevant sections in the documentation
- [ ] Added unit tests for new algorithms, conversion methods, etc.
- [ ] Added reference image tests for new plotting functions, recipes, visual options, etc.
# MakieCore

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://JuliaPlots.github.io/MakieCore.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://JuliaPlots.github.io/MakieCore.jl/dev)
[![Build Status](https://github.com/JuliaPlots/MakieCore.jl/workflows/CI/badge.svg)](https://github.com/JuliaPlots/MakieCore.jl/actions)
[![Coverage](https://codecov.io/gh/JuliaPlots/MakieCore.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/JuliaPlots/MakieCore.jl)
WGLMakie is a WebGL backend for the [Makie.jl](https://www.github.com/JuliaPlots/Makie.jl) plotting package, implemented using Three.js.

Read the docs for Makie and it's backends [here](http://makie.juliaplots.org/dev)

# Usage

Now, it should just work like Makie:

```julia
using WGLMakie
scatter(rand(4))
```

In the REPL, this will open a browser tab, that will refresh on a new display.
In VSCode, this should open in the plotpane.
You can also embed plots in a JSServe webpage:

```julia
function dom_handler(session, request)
    return DOM.div(
        DOM.h1("Some Makie Plots:"),
        meshscatter(1:4, color=1:4),
        meshscatter(1:4, color=rand(RGBAf, 4)),
        meshscatter(1:4, color=rand(RGBf, 4)),
        meshscatter(1:4, color=:red),
        meshscatter(rand(Point3f, 10), color=rand(RGBf, 10)),
        meshscatter(rand(Point3f, 10), marker=Pyramid(Point3f(0), 1f0, 1f0)),
    )
end
isdefined(Main, :app) && close(app)
app = JSServe.Server(dom_handler, "127.0.0.1", 8082)
```

## Sponsors

<img src="https://github.com/JuliaPlots/Makie.jl/blob/master/assets/BMBF_gefoerdert_2017_en.jpg?raw=true" width="300"/>
Förderkennzeichen: 01IS10S27, 2020
The OpenGL backend for [Makie](https://github.com/JuliaPlots/Makie.jl)

Read the docs for Makie and its backends [here](http://makie.juliaplots.org/dev)

## Issues
Please file all issues in [Makie.jl](https://github.com/JuliaPlots/Makie.jl/issues/new), and mention GLMakie in the issue text!


## Troubleshooting OpenGL

If you get any error loading GLMakie, it likely means, you don't have an OpenGL capable Graphic Card, or you don't have an OpenGL 3.3 capable video driver installed.
Note, that most GPUs, even 8 year old integrated ones, support OpenGL 3.3.

On Linux, you can find out your OpenGL version with:
`glxinfo | grep "OpenGL version"`

If you're using an AMD or Intel gpu on linux, you may run into [GLFW#198](https://github.com/JuliaGL/GLFW.jl/issues/198).

If you're on a headless server, you still need to install x-server and
proper graphics drivers.

You can find instructions to set that up in:

https://nextjournal.com/sdanisch/GLMakie-nogpu
And for a headless github action:

https://github.com/JuliaPlots/Makie.jl/blob/master/.github/workflows/glmakie.yaml
If none of these work for you, there is also a Cairo and WebGL backend
for Makie which you can use:

https://github.com/JuliaPlots/Makie.jl/tree/master/CairoMakie.

https://github.com/JuliaPlots/Makie.jl/tree/master/WGLMakie.

If you get an error pointing to [GLFW.jl](https://github.com/JuliaGL/GLFW.jl), please look into the existing [GLFW issues](https://github.com/JuliaGL/GLFW.jl/issues), and also google for those errors. This is then very likely something that needs fixing in the  [glfw c library](https://github.com/glfw/glfw) or in the GPU drivers.


## WSL setup or X-forwarding

From: https://github.com/Microsoft/WSL/issues/2855#issuecomment-358861903

WSL runs OpenGL alright, but it is not a supported scenario.
From a clean Ubuntu install from the store do:

```
sudo apt install ubuntu-desktop mesa-utils
export DISPLAY=localhost:0
glxgears
```

On the Windows side:

1) install [VcXsrv](https://sourceforge.net/projects/vcxsrv/)
2) choose multiple windows -> display 0 -> start no client -> disable native opengl

Troubleshooting:

1.)  install: `sudo apt-get install -y xorg-dev mesa-utils xvfb libgl1 freeglut3-dev libxrandr-dev libxinerama-dev libxcursor-dev libxi-dev libxext-dev`

2.) WSL has some problems with passing through localhost, so one may need to use: `export DISPLAY=192.168.178.31:0`, with the local ip of the pcs network adapter, which runs VcXsrv

3.) One may need `mv /opt/julia-1.5.2/lib/julia/libstdc++.so.6 /opt/julia-1.5.2/lib/julia/libcpp.backup`, another form of [GLFW#198](https://github.com/JuliaGL/GLFW.jl/issues/198)
# CairoMakie

The Cairo Backend for Makie

Read the docs for Makie and it's backends [here](http://makie.juliaplots.org/stable)


## Issues

Please file all issues in [Makie.jl](https://github.com/JuliaPlots/Makie.jl/issues/new), and mention CairoMakie in the issue text.

## Limitations

CairoMakie is intended as a backend for static vector graphics at publication quality. Therefore, it does not support the interactive features of GLMakie and is slower when visualizing large amounts data. 3D plots are currently not available because of the inherent limitations of 2D vector graphics.

## Saving

Makie overloads the FileIO interface, so you can save a Scene `scene` as `save("filename.extension", scene)`. CairoMakie supports saving to PNG, PDF, SVG and EPS.

You can scale the size of the output figure, without changing its appearance by passing keyword arguments to `save`. PNGs can be scaled by `px_per_unit` (default 1) and vector graphics (SVG, PDF, EPS) can be scaled by `pt_per_unit`.

```julia
save("plot.svg", scene, pt_per_unit = 0.5) # halve the dimensions of the resulting SVG
save("plot.png", scene, px_per_unit = 2) # double the resolution of the resulting PNG
```

## Using CairoMakie with Gtk.jl

You can render onto a GtkCanvas using Gtk, and use that as a display for your scenes.

```julia
using Gtk, CairoMakie

canvas = @GtkCanvas()
window = GtkWindow(canvas, "Makie", 500, 500)

function drawonto(canvas, figure)
    @guarded draw(canvas) do _
        scene = figure.scene
        resize!(scene, Gtk.width(canvas), Gtk.height(canvas))
        screen = CairoMakie.CairoScreen(scene, Gtk.cairo_surface(canvas), getgc(canvas), nothing)
        CairoMakie.cairo_draw(screen, scene)
    end
end

fig, ax, pl = heatmap(rand(50, 50)) # or something

drawonto(canvas, fig)
show(canvas); # trigger rendering
```
@def title = "Home"
@def order = 0
@def frontpage = true

# Welcome to Makie!

Makie is a data visualization ecosystem for the [Julia](https://julialang.org/) programming language, with high performance and extensibility.
It is available for Windows, Mac and Linux.

## Example

~~~
<input id="hidecode" class="hidecode" type="checkbox">
~~~
```julia:lorenz
using GLMakie
GLMakie.activate!() # hide

Base.@kwdef mutable struct Lorenz
    dt::Float64 = 0.01
    σ::Float64 = 10
    ρ::Float64 = 28
    β::Float64 = 8/3
    x::Float64 = 1
    y::Float64 = 1
    z::Float64 = 1
end

function step!(l::Lorenz)
    dx = l.σ * (l.y - l.x)
    dy = l.x * (l.ρ - l.z) - l.y
    dz = l.x * l.y - l.β * l.z
    l.x += l.dt * dx
    l.y += l.dt * dy
    l.z += l.dt * dz
    Point3f(l.x, l.y, l.z)
end

attractor = Lorenz()

points = Observable(Point3f[])
colors = Observable(Int[])

set_theme!(theme_black())

fig, ax, l = lines(points, color = colors,
    colormap = :inferno, transparency = true,
    axis = (; type = Axis3, protrusions = (0, 0, 0, 0),
        viewmode = :fit, limits = (-30, 30, -30, 30, 0, 50)))

record(fig, "lorenz.mp4", 1:120) do frame
    for i in 1:50
        push!(points[], step!(attractor))
        push!(colors[], frame)
    end
    ax.azimuth[] = 1.7pi + 0.3 * sin(2pi * frame / 120)
    notify.((points, colors))
    l.colorrange = (0, frame)
end
set_theme!() # hide
```
~~~
<label for="hidecode" class="hidecode"></label>
~~~

\video{lorenz, autoplay = true}

## Installation and Import

Add one or more of the Makie backend packages [`GLMakie.jl`](/documentation/backends/glmakie/) (OpenGL), [`CairoMakie.jl`](/documentation/backends/cairomakie/) (Cairo), or [`WGLMakie.jl`](/documentation/backends/wglmakie/) (WebGL), [`RPRMakie`](/documentation/backends/rprmakie/) (RadeonProRender) using Julia's inbuilt package manager. Each backend re-exports `Makie` so there's no need to install it separately.

```julia
]add GLMakie
using GLMakie
```

To switch to a different backend, for example `CairoMakie`, call `CairoMakie.activate!()`.

## First Steps

@@box-container
  @@box
    ~~~<a class="boxlink" href="tutorials/basic-tutorial/">~~~
    @@title Basic Tutorial @@
    @@box-content
      @@description
      Learn the basics of plotting with Makie.
      @@
      ~~~
      <img src="/assets/basic_tutorial_example.png">
      ~~~
    @@
    ~~~</a>~~~
  @@

  @@box
    ~~~<a class="boxlink" href="tutorials/layout-tutorial/">~~~
    @@title Layout Tutorial @@
    @@box-content
      @@description
      Check out how to make complex plots and layouts.
      @@
      ~~~
      <img src="/assets/tutorials/layout-tutorial/code/output/final_result.png">
      ~~~
    @@
    ~~~</a>~~~
  @@

  @@box
    ~~~<a class="boxlink" href="examples/plotting_functions/">~~~
    @@title Plot Examples @@
    @@box-content
      @@description
      Have a look at this list of examples for the available plotting functions.
      @@
      ~~~
      <img src="/assets/examples/plotting_functions/heatmap/code/output/mandelbrot_heatmap.png">
      ~~~
    @@
    ~~~</a>~~~
  @@

@@

## Makie Ecosystem

There are four backends, each of which has particular strengths. You can switch between backends at any time.

@@box-container
  @@box
    ~~~<a class="boxlink" href="/documentation/backends/glmakie/">~~~
    @@title GLMakie.jl@@
    @@box-content
      @@description
      GPU-powered, interactive 2D and 3D plotting in standalone `GLFW.jl` windows.
      @@
      ~~~
      <img src="/assets/surface_example.png">
      ~~~
    @@
    ~~~</a>~~~
  @@

  @@box
    ~~~<a class="boxlink" href="/documentation/backends/cairomakie/">~~~
    @@title CairoMakie.jl @@
    @@box-content
      @@description
      `Cairo.jl` based, non-interactive 2D backend for publication-quality vector graphics.
      @@
      ~~~
      <img src="/assets/density_example.png">
      ~~~
    @@
    ~~~</a>~~~
  @@

  @@box
    ~~~<a class="boxlink" href="/documentation/backends/wglmakie/">~~~
    @@title WGLMakie.jl @@
    @@box-content
      @@description
      WebGL-based interactive 2D and 3D plotting that runs within browsers.
      @@
      ~~~
      <img src="/assets/wireframe_example.png">
      ~~~
    @@
    ~~~</a>~~~
  @@
  @@box
    ~~~<a class="boxlink" href="documentation/backends/rprmakie/">~~~
    @@title RPRMakie.jl @@
    @@box-content
      @@description
      Backend using RadeonProRender for raytracing Makie scenes.
      @@
      ~~~
      <img src="/assets/topographie.png">
      ~~~
    @@
    ~~~</a>~~~
  @@
@@

The differences between backends are explained in more details under \myreflink{Backends}.

### Extensions and Resources

These packages and sites are maintained by third parties. If you install packages, keep an eye on version conflicts or downgrades as the Makie ecosystem is developing quickly so things break occasionally.

@@box-container
  @@box
    ~~~<a class="boxlink" href="https://github.com/JuliaPlots/AlgebraOfGraphics.jl/">~~~
    @@title AlgebraOfGraphics.jl @@
    @@box-content
      @@description
      Grammar-of-graphics style plotting, inspired by ggplot2.
      @@
      ~~~
      <img src="/assets/aog_example.png">
      ~~~
    @@
    ~~~</a>~~~
  @@

  @@box
    ~~~<a class="boxlink" href="https://lazarusa.github.io/BeautifulMakie/">~~~
    @@title Beautiful Makie @@
    @@box-content
      @@description
      This third-party gallery contains many advanced examples.
      @@
      ~~~
      <img src="/assets/beautifulmakie_example.png">
      ~~~
    @@
    ~~~</a>~~~
  @@

  @@box
    ~~~<a class="boxlink" href="https://github.com/JuliaPlots/GraphMakie.jl">~~~
    @@title GraphMakie.jl @@
    @@box-content
      @@description
      Graphs with two- and three-dimensional layout algorithms.
      @@
      ~~~
      <img src="/assets/graphmakie.png">
      ~~~
    @@
    ~~~</a>~~~
  @@

  @@box
    ~~~<a class="boxlink" href="https://github.com/JuliaPlots/GeoMakie.jl">~~~
    @@title GeoMakie.jl @@
    @@box-content
      @@description
      Geographic plotting utilities including projections.
      @@
      ~~~
      <img src="/assets/geomakie_example.png">
      ~~~
    @@
    ~~~</a>~~~
  @@
@@


## Citing Makie

If you use Makie for a scientific publication, please cite [our JOSS paper](https://joss.theoj.org/papers/10.21105/joss.03349) the following way:

> Danisch & Krumbiegel, (2021). Makie.jl: Flexible high-performance data visualization for Julia. Journal of Open Source Software, 6(65), 3349, https://doi.org/10.21105/joss.03349

You can use the following BibTeX entry:

```
@article{DanischKrumbiegel2021,
  doi = {10.21105/joss.03349},
  url = {https://doi.org/10.21105/joss.03349},
  year = {2021},
  publisher = {The Open Journal},
  volume = {6},
  number = {65},
  pages = {3349},
  author = {Simon Danisch and Julius Krumbiegel},
  title = {Makie.jl: Flexible high-performance data visualization for Julia},
  journal = {Journal of Open Source Software}
}
```

## Getting Help

1. Use the REPL `?` help mode.
1. Click this link to open a preformatted topic on the [Julia Discourse Page](https://discourse.julialang.org/new-topic?title=Makie%20-%20Your%20question%20here&category=domain/viz&tags=Makie&body=You%20can%20write%20your%20question%20in%20this%20space.%0A%0ABefore%20asking%2C%20please%20take%20a%20minute%20to%20make%20sure%20that%20you%20have%20installed%20the%20latest%20available%20versions%20and%20have%20looked%20at%20%5Bthe%20most%20recent%20documentation%5D(http%3A%2Fmakie.juliaplots.org%2Fstable%2F)%20%3Ainnocent%3A). If you do this manually, please use the category Domain/Visualization and tag questions with `Makie` to increase their visibility.
1. For casual conversation about Makie and its development, have a look at the `#makie` channel in the [Julia Slack group](https://julialang.org/slack/). Please direct your usage questions to [Discourse](https://discourse.julialang.org/new-topic?title=Makie%20-%20Your%20question%20here&category=domain/viz&tags=Makie&body=You%20can%20write%20your%20question%20in%20this%20space.%0A%0ABefore%20asking%2C%20please%20take%20a%20minute%20to%20make%20sure%20that%20you%20have%20installed%20the%20latest%20available%20versions%20and%20have%20looked%20at%20%5Bthe%20most%20recent%20documentation%5D(http%3A%2Fmakie.juliaplots.org%2Fstable%2F)%20%3Ainnocent%3A) and not to Slack, to make questions and answers accessible to everybody.
1. For technical issues and bug reports, open an issue in the [Makie.jl](https://github.com/JuliaPlots/Makie.jl) repository which serves as the central hub for Makie and backend issues.
@def order = 3

# Documentation

{{list_folder documentation}}<!--
Add here global page variables to use throughout your website.
-->
+++
author = "Makie.jl"
mintoclevel = 2
frontpage = false
auto_code_path = true

# Add here files or directories that should be ignored by Franklin, otherwise
# these files might be copied and, if markdown, processed by Franklin which
# you might not want. Indicate directories by ending the name with a `/`.
# Base files such as LICENSE.md and README.md are ignored by default.
ignore = ["node_modules/", "misc/"]

# RSS (the website_{title, descr, url} must be defined to get RSS)
generate_rss = true
website_title = "Franklin Template"
website_descr = "Example website using Franklin"
prepath = get(ENV, "PREVIEW_FRANKLIN_PREPATH", "")
website_url = get(ENV, "PREVIEW_FRANKLIN_WEBSITE_URL", "makie.juliaplots.org")
+++

<!--
Add here global latex commands to use throughout your pages.
-->
\newcommand{\R}{\mathbb R}
\newcommand{\scal}[1]{\langle #1 \rangle}

<!-- myreflink{Basic Tutorial} expands to [Basic Tutorial](link_to_that) -->
\newcommand{\myreflink}[1]{[!#1](\reflink{!#1})}

\newcommand{\apilink}[1]{[`!#1`](/documentation/api_reference/#!#1)}@def order = 2

# Tutorials

@@box-container
  @@box
    ~~~<a class="boxlink" href="basic-tutorial/">~~~
    @@title Basic Tutorial @@
    @@box-content
      @@description
      Learn the basics of plotting with Makie.
      @@
      ~~~
      <img src="/assets/basic_tutorial_example.png">
      ~~~
    @@
    ~~~</a>~~~
  @@

  @@box
    ~~~<a class="boxlink" href="layout-tutorial/">~~~
    @@title Layout Tutorial @@
    @@box-content
      @@description
      Check out how to make complex plots and layouts.
      @@
      ~~~
      <img src="/assets/tutorials/layout-tutorial/code/output/final_result.png">
      ~~~
    @@
    ~~~</a>~~~
  @@

  @@box
    ~~~<a class="boxlink" href="aspect-tutorial/">~~~
    @@title Aspect Tutorial @@
    @@box-content
      @@description
      How to deal with fixed sizes and aspect ratios.
      @@
      ~~~
      <img src="/assets/tutorials/aspect-tutorial/code/output/aspect_tutorial_example.png">
      ~~~
    @@
    ~~~</a>~~~
  @@

  @@box
    ~~~<a class="boxlink" href="scenes/">~~~
    @@title Scene Tutorial @@
    @@box-content
      @@description
      Explains the Scene type, projections and the scene graph.
      @@
      ~~~
      <video src="/assets/lego_walk.mp4">
      ~~~
    @@
    ~~~</a>~~~
  @@
@@

## Videos

@@box-container
  @@box
    ~~~<a class="boxlink" href="https://www.youtube.com/watch?v=AtqXFfaMZqo">~~~
    @@title MakieLayout Deep Dive @@
    @@box-content
      @@description
      A tour of some of the more complex layout features of Makie.
      @@
      ~~~
      <div class="youtube-container">
      <iframe width="560" height="315" src="https://www.youtube.com/embed/AtqXFfaMZqo?controls=0" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe>
      </div>
      ~~~
    @@
    ~~~</a>~~~
  @@

  @@box
    ~~~<a class="boxlink" href="https://www.youtube.com/watch?v=L-gyDvhjzGQ">~~~
    @@title Animations & Interaction @@
    @@box-content
      @@description
      How to create animations and interactive applications in Makie.
      @@
      ~~~
      <div class="youtube-container">
      <iframe width="560" height="315" src="https://www.youtube.com/embed/L-gyDvhjzGQ?controls=0" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe>
      </div>
      ~~~
    @@
    ~~~</a>~~~
  @@

@@
@def title = "Search"
@def hidden = true

## Search

Number of results found: ~~~<span id="resultCount"></span>~~~

~~~
<div id="searchResults"></div>
~~~@def title = "404"
@def hidden = true

~~~
<div style="margin-top: 40px; font-size: 40px; text-align: center;">

  <br>

  <div style="font-weight: bold;">
    404
  </div>

  <br>
  <br>

    The requested page was not found

  <br>
  <br>
  <br>
  <br>

  <div style="margin-bottom: 300px; font-size: 24px">
    <a href="/">Click here</a> to go back to the homepage.
  </div>

</div>
~~~
@def order = 1

# Examples

{{list_folder examples}}# Inspecting Data

Makie provides a data inspection tool via `DataInspector(x)` where x can be a
figure, axis or scene. With it you get a floating tooltip with relevant
information for various plots by hovering over one of its elements.

By default the inspector will be able to pick any plot other than `text` and
`volume` based plots. If you wish to ignore a plot, you can set its attribute
`plot.inspectable[] = false`. With that the next closest plot (in range) will be
picked.

## Attributes of `DataInspector`

The `inspector = DataInspector(fig)` contains the following attributes:

- `range = 10`: Controls the snapping range for selecting an element of a plot.
- `enabled = true`: Disables inspection of plots when set to false. Can also be adjusted with `enable!(inspector)` and `disable!(inspector)`.
- `text_padding = Vec4f(5, 5, 3, 3)`: Padding for the box drawn around the tooltip text. (left, right, bottom, top)
- `text_align = (:left, :bottom)`: Alignment of text within the tooltip. This does not affect the alignment of the tooltip relative to the cursor.
- `textcolor = :black`: Tooltip text color.
- `textsize = 20`: Tooltip text size.
- `font = "Dejavu Sans"`: Tooltip font.
- `background_color = :white`: Background color of the tooltip.
- `outline_color = :grey`: Outline color of the tooltip.
- `outline_linestyle = nothing`: Linestyle of the tooltip outline.
- `outline_linewidth = 2`: Linewidth of the tooltip outline.
- `indicator_color = :red`: Color of the selection indicator.
- `indicator_linewidth = 2`: Linewidth of the selection indicator.
- `indicator_linestyle = nothing`: Linestyle of the selection indicator
- `tooltip_align = (:center, :top)`: Default position of the tooltip relative to the cursor or current selection. The real align may adjust to keep the tooltip in view.
- `tooltip_offset = Vec2f(20)`: Offset from the indicator to the tooltip.
- `depth = 9e3`: Depth value of the tooltip. This should be high so that the tooltip is always in front.
- `priority = 100`: The priority of creating a tooltip on a mouse movement or scrolling event.

## Extending the `DataInspector`

The inspector implements tooltips for primitive plots and a few non-primitive
(i.e. a recipe) plots. All other plots will fall back to tooltips of their
hovered child. While this means that most plots have a tooltip it also means
many may not have a fitting one. If you wish to implement a more fitting tooltip
for non-primitive plot you may do so by creating a method

```julia
function show_data(inspector::DataInspector, my_plot::MyPlot, idx, primitive_child::SomePrimitive)
    ...
end
```

Here `my_plot` is the plot you want to create a custom tooltip for,
`primitive_child` is one of the primitives your plot is made from and `idx` is
the index into that primitive plot. The latter two are the result from
`pick_sorted` at the mouseposition. In general you will need to adjust `idx` to
be useful for `MyPlot`.

Let's take a look at the `BarPlot` overload, which also powers `hist`. It
contains two primitive plots - `Mesh` and `Lines`. The `idx` from picking a
`Mesh` is based on vertices, which there are four per rectangle. From `Lines` we
get an index based on the end point of the line. To draw the outline of a
rectangle we need 5 points and a seperator, totaling 6. We thus implement

```julia
import Makie: show_data

function show_data(inspector::DataInspector, plot::BarPlot, idx, ::Lines)
    return show_barplot(inspector, plot, div(idx-1, 6)+1)
end

function show_data(inspector::DataInspector, plot::BarPlot, idx, ::Mesh)
    return show_barplot(inspector, plot, div(idx-1, 4)+1)
end
```

to map the primitive `idx` to one relevant for `BarPlot`. With this we can now
get the position of the hovered bar with `plot[1][][idx]`. To align the tooltip
to the selection we need to compute the relevant position in screen space and
update the tooltip position.

```julia
using Makie: parent_scene, shift_project, update_tooltip_alignment!, position2string

function show_barplot(inspector::DataInspector, plot::BarPlot, idx)
    # All the attributes of DataInspector are here
    a = inspector.plot.attributes

    # Get the scene BarPlot lives in
    scene = parent_scene(plot)

    # Get the hovered world-space position
    pos = plot[1][][idx]
    # project to screen space and shift it to be correct on the root scene
    proj_pos = shift_project(scene, to_ndim(Point3f, pos, 0))
    # anchor the tooltip at the projected position
    update_tooltip_alignment!(inspector, proj_pos)

    # Update the final text of the tooltip
    # position2string is just an `@sprintf`
    a._display_text[] = position2string(pos)
    # Show the tooltip
    a._visible[] = true

    # return true to indicate that we have updated the tooltip
    return true
end
```

Next we need to mark the rectangle we are hovering. In this case we can use the
rectangles which `BarPlot` passes to `Poly`, i.e. `plot.plots[1][1][][idx]`. The
`DataInspector` contains some functionality for keeping track of temporary plots,
so we can plot the indicator to the same `scene` that `BarPlot` uses. Doing so
results in

```julia
using Makie:
    parent_scene, shift_project, update_tooltip_alignment!, position2string,
    clear_temporary_plots!


function show_barplot(inspector::DataInspector, plot::BarPlot, idx)
    a = inspector.plot.attributes
    scene = parent_scene(plot)

    pos = plot[1][][idx]
    proj_pos = shift_project(scene, to_ndim(Point3f, pos, 0))
    update_tooltip_alignment!(inspector, proj_pos)

    # Get the rectangle BarPlot generated for Poly
    # `_bbox2D` is a observable meant for saving a `Rect2` indicator. There is also
    # a `_bbox3D`. Here we keep `_bbox2D` updated and use it as a source for
    # our custom indicator.
    a._bbox2D[] = plot.plots[1][1][][idx]
    a._model[] = plot.model[]

    # Only plot the indicator once. It'll be updated via `_bbox2D`.
    if inspector.selection != plot
        # Clear any old temporary plots (i.e. other indicators like this)
        # this also updates inspector.selection.
        clear_temporary_plots!(inspector, plot)

        # create the indicator using a bunch of the DataInspector attributes.
        # Note that temporary plots only cleared when a new one is created. To
        # control whether indicator is visible or not `a._bbox_visible` is set
        # instead, so it should be in any custom indicator like this.
        p = wireframe!(
            scene, a._bbox2D, model = a._model, color = a.indicator_color,
            strokewidth = a.indicator_linewidth, linestyle = a.indicator_linestyle,
            visible = a._bbox_visible, inspectable = false
        )

        # Make sure this draws on top
        translate!(p, Vec3f(0, 0, a.depth[]))

        # register this indicator for later cleanup.
        push!(inspector.temp_plots, p)
    end

    a._display_text[] = position2string(pos)
    a._visible[] = true

    # Show our custom indicator
    a._bbox_visible[] = true
    # Don't show the default screen space indicator
    a._px_bbox_visible[] = false

    return true
end
```

which finishes the implementation of a custom tooltip for `BarPlot`.
# Events

Interactive backends such as `GLMakie` and `WGLMakie` pass events to an `Events`
struct in the currently active scene. There they trigger a number of slightly
modified observables of type `PriorityObservable`. By reacting to those one can
build up custom interactions.

## PriorityObservables

Much like the name suggests a `PriorityObservable` adds a priority to its
listeners. Furthermore it allows for each listener to stop execution of lower
priority listeners by returning `Consume(true)` or simply `Consume()`. Every
other return value will be handled as `Consume(false)` meaning that the 
listener does not block other listeners.

To understand how a `PriorityObserable` works you may try this example:

```julia
using Makie: PriorityObservable

po = PriorityObservable(0)

on(po, priority = -1) do x
    println("Low priority: $x")
end
po[] = 1

on(po, priority = 0) do x
    println("Medium blocking priority: $x")
    return Consume()
end
po[] = 2

on(po, priority = 1) do x
    println("High Priority: $x")
    return Consume(false)
end
po[] = 3
```

With only the first listener connected you should see `Low priority: 1` getting
printed. With the second you should only see a reaction from the medium priority
listener, as it blocks execution of the lower priority by returning `Consume()`. 
With all three connected you should see two prints, first the high priority, 
than the medium priority.

## The Events struct

Events from the backend are stored in `PriorityObservables` within the `Events`
struct. You can access it from any scene via `events(scene)` and via
`events(figure.scene)` or `events(axis.scene)` from a figure or axis. The struct
contains the following fields you may react to.

- `window_area::PriorityObservable{Rect2i}`: Contains the current size of the window in pixels.
- `window_dpi::PriorityObservable{Float64}`: Contains the DPI of the window.
- `window_open::PriorityObservable{Bool}`: Contains `true` as long as the window is open.
- `hasfocus::PriorityObservable{Bool}`: Contains `true` if the window is focused (in the foreground).
- `entered_window::PriorityObservable{Bool}`: Contains true if the mouse is within the window (regarless of whether it is focused), i.e. when it is hovered.
- `mousebutton::PriorityObservable{MouseButtonEvent}`: Contains the most recent `MouseButtonEvent` which holds the relevant `button::Mouse.Button` and `action::Mouse.Action`.
- `mousebuttonstate::Set{Mouse.Button}`: Contains all currently pressed mouse buttons.
- `mouseposition::PriorityObservable{NTuple{2, Float64}}`: Contains the most recent cursor position in pixel units relative to the root scene/window.
- `scroll::PriorityObservable{NTuple{2, Float64}}`: Contains the most recent scroll offset.
- `keyboardbutton::PriorityObservable{KeyEvent}`: Contains the most recent `KeyEvent` which holds the relevant `key::Keyboard.Button` and `action::Keyboard.Action`.
- `keyboardstate::PriorityObservable{Keyboard.Button}`: Contains all currently pressed keys.
- `unicode_input::PriorityObservable{Char}`: Contains the most recently typed character.
- `dropped_files::PriorityObservable{Vector{String}}`: Contains a list of filepaths to a collection files dragged into the window.

## Mouse Interaction

There are three mouse events one can react to:

- `events.mousebutton` which holds a `MouseButtonEvent` with relevant `button` and `action`
- `events.mouseposition` which holds the current cursor position relative to the window as `NTuple{2, Float64}` in pixel
- `events.scroll` which holds an `NTuple{2, Float64}` of the last scroll change

For example, we can react to a left click with:

```julia
on(events(fig).mousebutton, priority = 0) do event
    if event.button == Mouse.left
        if event.action == Mouse.press
            # do something
        else
            # do something else when the mouse button is released
        end
    end
    # Do not consume the event
    return Consume(false)
end
```

There are also a bunch of convenience function for mouse interactions:

- `hovered_scene()` returns the currently hovered scene.
- `mouseposition_px([scene])` returns the cursor position in pixel units relative to the given or currently hovered scene.
- `mouseposition([scene])` returns the cursor position of the given or hovered scene in data coordinates.
- `mouseover(scene, plots...)` returns true if the mouse is over any of the given plots.
- `mouse_selection(scene[, range])` returns the plot under your cursor.
- `onpick(f, scene, plots...[; range=1])` allows you to define a callback `f` when one of the given plots is hovered (or the closest within range of the cursor).
- `pick(scene, position[, range])` returns the plot under the given position (or the closest within range)
- `Makie.pick_closest(scene, position, range)` returns all plots within range, sorted by distance to the given position.

## Keyboard Interaction

You can use `events.keyboardbutton` to react to a `KeyEvent` and
`events.unicode_input` to react to specific characters being typed.

For example we may react to specific keys being pressed or held with

```julia
on(events(fig).keyboardbutton) do event
    if event.action in (Keyboard.press, Keyboard.repeat)
        event.key == Keyboard.left   && move_left()
        event.key == Keyboard.up     && move_up()
        event.key == Keyboard.right  && move_right()
        event.key == Keyboard.down   && move_down()
    end
    # Let the event reach other listeners
    return Consume(false)
end
```

## The `ispressed` function

`ispressed(scene, x)` can be used to check if a button, collection of buttons or
a boolean expression of buttons is pressed.

```julia
# Check if the "a" key is pressed
ispressed(scene, Keyboard.a)

# Check if the left mouse button is pressed
ispressed(scene, Mouse.left)

# check if the left mouse button, left shift and left control are pressed.
# This also works with Vectors and Sets.
ispressed(scene, (Mouse.left, Keyboard.left_shift, Keyboard.left_control))

# check if the left mouse button and either left or right control are pressed
ispressed(scene, Mouse.left & (Keyboard.left_control | Keyboard.right_control))

# check if the "c" key is pressed without either shift key
ispressed(scene, Keyboard.c & !(Keyboard.left_shift | Keyboard.right_shift))
```

In all the examples above it is irrelevant how many buttons are pressed as long 
as the specified subset is pressed (or the boolean expression matches). You can
also check for exact matches by wrapping a button, collection or expression in
`Exclusively`. For example

```julia
# check if exclusively "a" and left shift are pressed (i.e. no other key or mouse button)
ispressed(scene, Exclusively(Keyboard.a & Keyboard.left_shift))
```

You can also pass true (false) to `ispressed`. In this case true (false) will
always be returned. This is useful for removing interactivity, i.e. always or 
never triggering something.

## Interactive Widgets

Makie has a couple of useful interactive widgets like sliders, buttons and menus, which you can read about in the \myreflink{Layoutables} section.

## Recording Animations with Interactions

You can record a `Scene` while you're interacting with it.
Just use the \apilink{record} function (also see the \myreflink{Animations} page) and allow interaction by `sleep`ing in the loop.

In this example, we sample from the Scene `scene` for 10 seconds, at a rate of 10 frames per second.

```julia
fps = 10
record(scene, "test.mp4"; framerate = fps) do io
    for i = 1:100
        sleep(1/fps)
        recordframe!(io)
    end
end
```
# Headless

Makie can be used on headless systems (such as CI servers).
This page describes what is required to get different back ends working in headless systems.

## CairoMakie

For CairoMakie, there shouldn't be any difference in using it on a remote or locally.

## GLMakie

For [GLMakie](https://github.com/JuliaPlots/GLMakie.jl) you can either use X11 forwarding to render on the local
host or use [VirtualGL](https://www.virtualgl.org/) to render on the remote server.

### GLMakie with X11 forwarding

In this scenario you need an X server on the remote and you will have to connect to the remote server with
```
ssh -X user@host
```
See [here](https://unix.stackexchange.com/questions/12755/how-to-forward-x-over-ssh-to-run-graphics-applications-remotely)
about more details about X11 forwarding.

### GLMakie with VirtualGL

The first step is to [install VirtualGL](https://cdn.rawgit.com/VirtualGL/virtualgl/2.6.3/doc/index.html#hd005) on the remote
server ([Linux only](https://virtualgl.org/Documentation/OSSupport)) and on the local client.
If you need to establish the connection to the server via a secondary intermediate server,
VirtualGL also needs to be installed there.
On the remote server you will need to [configure the VirtualGL server](https://cdn.rawgit.com/VirtualGL/virtualgl/2.6.5/doc/index.html#hd006).
Be sure to [check that the configuration is ok](https://cdn.rawgit.com/VirtualGL/virtualgl/2.6.5/doc/index.html#hd006002001).

After everything is set up, you can connect to the remote server via
```
/opt/VirtualGL/bin/vglconnect -s user@server
```
and then you will have to start julia via VirtualGL
```
/opt/VirtualGL/bin/vglrun julia
```

### GLMakie in CI

You can also use GLMakie on CI or servers without a GPU by using `xvfb` for software rendering.
This procedure is used in the [GLMakie tests](https://github.com/JuliaPlots/GLMakie.jl/blob/8954fc34354a09ceb11159a8e8e35429c05a710f/.github/workflows/ci.yml#L41-L42).

## WGLMakie

For WGLMakie, you can setup a server with JSServe and serve the content from a remote server.
This also works for creating interactive plots with Documenter.
Check out the [docs](http://juliaplots.org/WGLMakie.jl/stable/) for more details about this.

If you want to use WGLMakie in VS Code on a remote server, you will have to forward the port
used by WGLMakie in order for the plot pane integration to work.
If you don't need to change the port on which WGLMakie,
you will just have to [forward](https://code.visualstudio.com/docs/remote/ssh#_forwarding-a-port-creating-ssh-tunnel) the 9284 port.

If you want to change the port on which WGLMakie runs on the remote, say `8081`, you will have to use the following
```julia
using JSServe

JSServe.configure_server!(listen_port=8081)
```
before any plotting commands with WGLMakie.

If you also need to use a different port than `8081` on the _local_ machine, say `8080`,
you will also need to set the `forwarded_port` like this:
```julia
using JSServe

JSServe.configure_server!(listen_port=8081, forwarded_port=8080)
```
# Theming

Makie allows you to change almost every visual aspect of your plots via attributes.
You can set attributes whenever you create an object, or you define a general style that is then used as the default by all following objects.

There are three functions you can use for that purpose:

```julia
set_theme!
update_theme!
with_theme
```

## set_theme!

You can call `set_theme!(theme; kwargs...)` to change the current default theme to `theme` and override or add attributes given by `kwargs`.
You can also reset your changes by calling `set_theme!()` without arguments.

Let's create a plot with the default theme:

\begin{examplefigure}{}
```julia
using CairoMakie
CairoMakie.activate!() # hide
Makie.inline!(true) # hide

function example_plot()
    f = Figure()
    for i in 1:2, j in 1:2
        lines(f[i, j], cumsum(randn(50)))
    end
    f
end

example_plot()
```
\end{examplefigure}

Now we define a theme which changes the default fontsize, activate it, and plot.

\begin{examplefigure}{}
```julia
fontsize_theme = Theme(fontsize = 10)
set_theme!(fontsize_theme)

example_plot()
```
\end{examplefigure}

This theme will be active until we call `set_theme!()`.

```julia:set_theme
set_theme!()
```

## update_theme!

If you have activated a theme already and want to update it partially, without removing the attributes not in the new theme, you can use `update_theme!`.

For example, you can first call `set_theme!(my_theme)` and later update font and fontsize with `update_theme!(font = "Arial", fontsize = 18)`, leaving all other settings intact.


## with_theme

Because it can be tedious to remember to switch themes off which you need only temporarily, there's the function `with_theme(f, theme)` which handles the resetting for you automatically, even if you encounter an error while running `f`.

\begin{examplefigure}{}
```julia
with_theme(fontsize_theme) do
    example_plot()
end
```
\end{examplefigure}

You can also pass additional keywords to add or override attributes in your theme:

\begin{examplefigure}{}
```julia
with_theme(fontsize_theme, fontsize = 25) do
    example_plot()
end
```
\end{examplefigure}

## Theming plot objects

You can theme plot objects by using their uppercase type names as a key in your theme.

\begin{examplefigure}{}
```julia
lines_theme = Theme(
    Lines = (
        linewidth = 4,
        linestyle = :dash,
    )
)

with_theme(example_plot, lines_theme)
```
\end{examplefigure}

## Theming layoutable objects

Every Layoutable such as `Axis`, `Legend`, `Colorbar`, etc. can be themed by using its type name as a key in your theme.

Here is how you could define a simple ggplot-like style for your axes:

\begin{examplefigure}{}
```julia
ggplot_theme = Theme(
    Axis = (
        backgroundcolor = :gray90,
        leftspinevisible = false,
        rightspinevisible = false,
        bottomspinevisible = false,
        topspinevisible = false,
        xgridcolor = :white,
        ygridcolor = :white,
    )
)

with_theme(example_plot, ggplot_theme)
```
\end{examplefigure}

## Cycles

Makie supports a variety of options for cycling plot attributes automatically.
For a plot object to use cycling, either its default theme or the currently active theme must have the `cycle` attribute set.

There are multiple ways to specify this attribute:

```julia
# You can either make a list of symbols
cycle = [:color, :marker]
# or map specific plot attributes to palette attributes
cycle = [:linecolor => :color, :marker]
# you can also map multiple attributes that should receive
# the same cycle attribute
cycle = [[:linecolor, :markercolor] => :color, :marker]
# nothing disables cycling
cycle = nothing # equivalent to cycle = []
```

### Covarying cycles

You can also construct a `Cycle` object directly, which additionally allows to set the `covary` keyword, that defaults to `false`. A cycler with `covary = true` cycles all attributes together, instead of cycling through all values of the first, then the second, etc.

```julia
# palettes: color = [:red, :blue, :green] marker = [:circle, :rect, :utriangle, :dtriangle]

cycle = [:color, :marker]
# 1: :red, :circle
# 2: :blue, :circle
# 3: :green, :circle
# 4: :red, :rect
# ...

cycle = Cycle([:color, :marker], covary = true)
# 1: :red, :circle
# 2: :blue, :rect
# 3: :green, :utriangle
# 4: :red, :dtriangle
# ...
```

### Manual cycling using `Cycled`

If you want to give a plot's attribute a specific value from the respective cycler, you can use the `Cycled` object.
The index `i` passed to `Cycled` is used directly to look up a value in the cycler that belongs to the attribute, and errors if no such cycler is defined.
For example, to access the third color in a cycler, instead of plotting three plots to advance the cycler, you can use `color = Cycled(3)`.

The cycler's internal counter is not advanced when using `Cycled` for any attribute, and only attributes with `Cycled` access the cycled values, all other usually cycled attributes fall back to their non-cycled defaults.

\begin{examplefigure}{}
```julia
using CairoMakie
CairoMakie.activate!() # hide
Makie.inline!(true) # hide

f = Figure()

Axis(f[1, 1])

# the normal cycle
lines!(0..10, x -> sin(x) - 1)
lines!(0..10, x -> sin(x) - 2)
lines!(0..10, x -> sin(x) - 3)

# manually specified colors
lines!(0..10, x -> sin(x) - 5, color = Cycled(3))
lines!(0..10, x -> sin(x) - 6, color = Cycled(2))
lines!(0..10, x -> sin(x) - 7, color = Cycled(1))

f
```
\end{examplefigure}

### Palettes

The attributes specified in the cycle are looked up in the axis' palette.
A single `:color` is both plot attribute as well as palette attribute, while `:color => :patchcolor` means that `plot.color` should be set to `palette.patchcolor`.
Here's an example that shows how density plots react to different palette options:

\begin{examplefigure}{}
```julia
using CairoMakie
CairoMakie.activate!() # hide
Makie.inline!(true) # hide

set_theme!() # hide

f = Figure(resolution = (800, 800))

Axis(f[1, 1], title = "Default cycle palette")

for i in 1:6
    density!(randn(50) .+ 2i)
end

Axis(f[2, 1],
    title = "Custom cycle palette",
    palette = (patchcolor = [:red, :green, :blue, :yellow, :orange, :pink],))

for i in 1:6
    density!(randn(50) .+ 2i)
end

set_theme!(Density = (cycle = [],))

Axis(f[3, 1], title = "No cycle")

for i in 1:6
    density!(randn(50) .+ 2i)
end

set_theme!() # hide

f
```
\end{examplefigure}

You can also theme global palettes via `set_theme!(palette = (color = my_colors, marker = my_markers))` for example.

## Special attributes

You can use the keys `rowgap` and `colgap` to change the default grid layout gaps.
# LaTeX

Makie can render LaTeX strings from the [LaTeXStrings.jl](https://github.com/stevengj/LaTeXStrings.jl) package using [MathTeXEngine.jl](https://github.com/Kolaru/MathTeXEngine.jl/).

This engine supports a subset of LaTeX's most used commands, which are rendered quickly enough for responsive use in GLMakie.

## Using L-strings

You can pass `LaTeXString` objects to almost any object with text labels. They are constructed using the `L` string macro prefix.
The whole string is interpreted as an equation if it doesn't contain an unescaped `$`.


```!
# hideall
using CairoMakie
```

\begin{examplefigure}{svg = true}
```julia
using CairoMakie
CairoMakie.activate!() # hide

f = Figure(fontsize = 18)

Axis(f[1, 1],
    title = L"\frac{x + y}{\sin(k^2)}",
    xlabel = L"\sum_a^b{xy}",
    ylabel = L"\sqrt{\frac{a}{b}}"
)

f
```
\end{examplefigure}

You can also mix math-mode and text-mode.
For [string interpolation](https://docs.julialang.org/en/v1/manual/strings/#string-interpolation) use `%$`instead of `$`:

\begin{examplefigure}{svg = true}
```julia
using CairoMakie
CairoMakie.activate!() # hide

f = Figure(fontsize = 18)
t = "text"
Axis(f[1,1], title=L"Some %$(t) and some math: $\frac{2\alpha+1}{y}$")

f
```
\end{examplefigure}


## Caveats

You can currently not change a text object which has been instantiated with a normal `String` input to use a `LaTeXString`. If you do this, the `LaTeXString` is converted to a `String` implicitly, losing its special properties. Instead it appears with `$` signs wrapped around.

Notice how the second axis loses its LaTeX title. You have to pass the L-string at axis construction.

\begin{examplefigure}{svg = true}
```julia
using CairoMakie
CairoMakie.activate!() # hide

f = Figure(fontsize = 18)

Axis(f[1, 1], title = L"\frac{x + y}{\sin(k^2)}")
ax2 = Axis(f[1, 2])
ax2.title = L"\frac{x + y}{\sin(k^2)}"

f
```
\end{examplefigure}


For implicitly created axes, you can do this via the `axis` keyword argument.


\begin{examplefigure}{svg = true}
```julia
scatter(randn(50, 2), axis = (; title = L"\frac{x + y}{\sin(k^2)}"))
```
\end{examplefigure}
# Backends

Makie is the frontend package that defines all plotting functions.
It is reexported by every backend, so you don't have to specifically install or import it.

There are four backends which concretely implement all abstract rendering capabilities defined in Makie:

| Package                                                        | Description                                                                           |
| :------------------------------------------------------------- | :------------------------------------------------------------------------------------ |
| [`GLMakie.jl`](/documentation/backends/glmakie/)       | GPU-powered, interactive 2D and 3D plotting in standalone `GLFW.jl` windows.          |
| [`CairoMakie.jl`](/documentation/backends/cairomakie/) | `Cairo.jl` based, non-interactive 2D (and some 3D) backend  for publication-quality vector graphics. |
| [`WGLMakie.jl`](/documentation/backends/wglmakie/)     | WebGL-based interactive 2D and 3D plotting that runs within browsers.                 |
| [`RPRMakie.jl`](/documentation/backends/rprmakie/)     | An experimental Ray tracing backend.                 |

### Activating Backends

You can activate any backend by `using` the appropriate package and calling its `activate!` function.

Example with WGLMakie:

```julia
using WGLMakie
WGLMakie.activate!()
```
# Layoutables

`Layoutables` are objects which can be added to a `Figure` or `Scene` and have their location and size controlled by a `GridLayout`. In of itself, a `Layoutable` is an abstract type.
A `Figure` has its own internal `GridLayout` and therefore offers simplified syntax for adding layoutables to it.
If you want to work with a bare `Scene`, you can attach a `GridLayout` to its pixel area.

!!! note
    A layout only controls an object's position or bounding box.
    A `Layoutable` can be controlled by the GridLayout of a Figure but not be added as a visual to the Figure.
    A `Layoutable` can also be added to a Scene without being inside any GridLayout, if you specify the bounding box yourself.

## Adding to a `Figure`

Here's one way to add a `Layoutable`, in this case an `Axis`, to a Figure.

\begin{examplefigure}{}
```julia
using CairoMakie
CairoMakie.activate!() # hide
Makie.inline!(true) # hide

f = Figure()
ax = Axis(f[1, 1])
f
```
\end{examplefigure}

## Specifying a boundingbox directly

Sometimes you just want to place a `Layoutable` in a specific location, without it being controlled by a dynamic layout.
You can do this by setting the `bbox` parameter, which is usually controlled by the layout, manually.
The boundingbox should be a 2D `Rect`, and can also be an Observable if you plan to change it dynamically.
The function `BBox` creates an `Rect2f`, but instead of passing origin and widths, you pass left, right, bottom and top boundaries directly.

Here's an example where two axes are placed manually:

\begin{examplefigure}{}
```julia
using CairoMakie
CairoMakie.activate!() # hide
Makie.inline!(true) # hide

f = Figure()
Axis(f, bbox = BBox(100, 300, 100, 500), title = "Axis 1")
Axis(f, bbox = BBox(400, 700, 200, 400), title = "Axis 2")
f
```
\end{examplefigure}
## Deleting layoutables

To remove layoutables from their layout and the figure or scene, use `delete!(layoutable)`.
# Animations

With Makie it is easy to create animated plots.
Animations work by making changes to data or plot attribute Observables and recording the changing figure frame by frame.
You can find out more about the Observables workflow on the \myreflink{Observables & Interaction} page.

## A simple example

To create an animation you need to use the \apilink{record} function.

First you create a `Figure`. Next, you pass a function that modifies this figure frame-by-frame to `record`.
Any changes you make to the figure or its plots will appear in the final animation.
You also need to pass an iterable which has as many elements as you want frames in your animation.
The function that you pass as the first argument is called with each element from this iterator
over the course of the animation.

As a start, here is how you can change the color of a line plot:

```julia:color_animation
using GLMakie
GLMakie.activate!() # hide
using Makie.Colors

fig, ax, lineplot = lines(0..10, sin; linewidth=10)

# animation settings
nframes = 30
framerate = 30
hue_iterator = range(0, 360, length=nframes)

record(fig, "color_animation.mp4", hue_iterator;
        framerate = framerate) do hue
    lineplot.color = HSV(hue, 1, 0.75)
end
nothing # hide
```

\video{color_animation}

Passing a function as the first argument is usually done with Julia's `do`-notation, which you might not be familiar with.
Instead of the above, we could also have written:

```julia
function change_function(hue)
    lineplot.color = HSV(hue, 1, 0.75)
end

record(change_function, fig, "color_animation.mp4", hue_iterator; framerate = framerate)
```

## File formats

Video files are created with [`FFMPEG.jl`](https://github.com/JuliaIO/FFMPEG.jl).
You can choose from the following file formats:

- `.mkv` (the default, doesn't need to convert)
- `.mp4` (good for web, widely supported)
- `.webm` (smallest file size)
- `.gif` (lowest quality with largest file size)

## Animations using `Observables`

Often, you want to animate a complex plot over time, and all the data that is displayed should be determined by the current time stamp.
Such a dependency is really easy to express with `Observables`.

We can save a lot of work if we create our data depending on a single time `Observable`, so we don't have to change every plot's data manually as the animation progresses.

Here is an example that plots two different functions.
The y-values of each depend on time and therefore we only have to change the time for both plots to change.
We use the convenient `@lift` macro which denotes that the `lift`ed expression depends on each Observable marked with a `$` sign.

```julia:time_animation
time = Observable(0.0)

xs = range(0, 7, length=40)

ys_1 = @lift(sin.(xs .- $time))
ys_2 = @lift(cos.(xs .- $time) .+ 3)

fig = lines(xs, ys_1, color = :blue, linewidth = 4,
    axis = (title = @lift("t = $(round($time, digits = 1))"),))
scatter!(xs, ys_2, color = :red, markersize = 15)

framerate = 30
timestamps = range(0, 2, step=1/framerate)

record(fig, "time_animation.mp4", timestamps;
        framerate = framerate) do t
    time[] = t
end
nothing # hide
```

\video{time_animation}

You can set most plot attributes equal to `Observable`s, so that you need only update
a single variable (like time) during your animation loop.

For example, to make a line with color dependent on time, you could write:

```julia:color_animation_2
time = Observable(0.0)
color_observable = @lift(RGBf($time, 0, 0))

fig = lines(0..10, sin, color = color_observable)

record(fig, "color_animation_2.mp4", timestamps; framerate = framerate) do t
    time[] = t
end
nothing # hide
```

\video{color_animation_2}

## Appending data with Observables

You can also append data to a plot during an animation.
Instead of passing `x` and `y` (or `z`) values separately,
it is better to make a `Observable` with a vector of `Point`s,
so that the number of `x` and `y` values can not go out of sync.

```julia:append_animation
points = Observable(Point2f[(0, 0)])

fig, ax = scatter(points)
limits!(ax, 0, 30, 0, 30)

frames = 1:30

record(fig, "append_animation.mp4", frames;
        framerate = 30) do frame
    new_point = Point2f(frame, frame)
    points[] = push!(points[], new_point)
end
nothing # hide
```

\video{append_animation}

## Animating a plot "live"

You can animate a live plot easily using a loop.
Update all `Observables` that you need and then add a short sleep interval so that the display can refresh:

```julia
points = Observable(Point2f[randn(2)])

fig, ax = scatter(points)
limits!(ax, -4, 4, -4, 4)

fps = 60
nframes = 120

for i = 1:nframes
    new_point = Point2f(randn(2))
    points[] = push!(points[], new_point)
    sleep(1/fps) # refreshes the display!
end
nothing # hide
```# Plot Method Signatures

Makie offers a simple but powerful set of methods for each plotting function, which allow you to easily create and manipulate the most common aspects of a figure.

Each plot object like `Scatter` has two plotting functions associated to it, a non-mutating version (`scatter`) and a mutating version (`scatter!`).
These functions have different methods that behave slightly differently depending on the first argument.

Here's a short list before we show each version in more detail.
We use `Scatter` as our example, but the principles apply to every plot type.

## Non-Mutating

The non-mutating methods create and return something in addition to the plot object, either a figure with an axis in default position, or an axis at a given GridPosition or GridSubposition.

```julia
scatter(args...; kwargs...) -> ::FigureAxisPlot
scatter(gridposition, args...; kwargs...) -> ::AxisPlot
```

`FigureAxisPlot` is just a collection of the new figure, axis and plot.
For convenience it has the same display overload as `Figure`, so that `scatter(args...)` displays a plot without further work.
It can be destructured at assignment like `fig, ax, plotobj = scatter(args...)`.

`AxisPlot` is a collection of a new axis and plot.
It has no special display overload but can also be destructured like `ax, plotobj = scatter(gridposition, args...)`.

### Special Keyword Arguments

Methods that create an `AxisPlot` accept a special-cased `axis` keyword, where you can pass a dict-like object containing keyword arguments that should be passed to the created axis.
Methods that create a `FigureAxisPlot` additionally accept a special cased `figure` keyword, where you can pass a dict-like object containing keyword arguments that should be passed to the created figure.

All other keyword arguments are passed as attributes to the plotting function.

Here are two examples with the scatter function (take care to create single-argument NamedTuples correctly, for example with a trailing comma):

\begin{examplefigure}{}
```julia
using GLMakie
GLMakie.activate!() # hide
# FigureAxisPlot takes figure and axis keywords
fig, ax, p = lines(cumsum(randn(1000)),
    figure = (resolution = (1000, 600),),
    axis = (ylabel = "Temperature",),
    color = :red)

# AxisPlot takes axis keyword
lines(fig[2, 1], cumsum(randn(1000)),
    axis = (xlabel = "Time (sec)", ylabel = "Stock Value"),
    color = :blue)

fig
```
\end{examplefigure}

## Mutating

The mutating methods always just return a plot object.
If no figure is passed, the `current_figure()` is used, if no axis or scene is given, the `current_axis()` is used.

```julia
scatter!(args...; kwargs...) -> ::Scatter
scatter!(figure, args...; kwargs...) -> ::Scatter
scatter!(gridposition, args...; kwargs...) -> ::Scatter
scatter!(axis, args...; kwargs...) -> ::Scatter
scatter!(scene, args...; kwargs...) -> ::Scatter
```

## GridPositions

In the background, each `Figure` has a `GridLayout` from [GridLayoutBase.jl](https://github.com/jkrumbiegel/GridLayoutBase.jl), which takes care of layouting plot elements nicely.
For convenience, you can index into a figure multiple times to refer to nested grid positions, which makes it easy to quickly assemble complex layouts.

For example, `fig[1, 2]` creates a `GridPosition` referring to row 1 and column 2, while `fig[1, 2][3, 1:2]` creates a `GridSubposition` that refers to row 3 and columns 1 to 2 in a nested GridLayout which is located at row 1 and column 2.
The link to the Figure is in the parent field of the top layout.

### With Non-Mutating Plotting Functions

Using the non-mutating plotting functions with GridPositions creates new axes at the given locations.
If a GridLayout along the nesting levels doesn't exist, yet, it is created automatically for convenience.

\begin{examplefigure}{}
```julia
using GLMakie
GLMakie.activate!() # hide
fig = Figure()

# first row, first column
scatter(fig[1, 1], 1.0..10, sin)

# first row, second column
lines(fig[1, 2], 1.0..10, sin)

# first row, third column, then nested first row, first column
lines(fig[1, 3][1, 1], cumsum(randn(1000)), color = :blue)

# first row, third column, then nested second row, first column
lines(fig[1, 3][2, 1], cumsum(randn(1000)), color = :red)

# second row, first to third column
ax, hm = heatmap(fig[2, 1:3], randn(30, 10))

# across all rows, new column after the last one
fig[:, end+1] = Colorbar(fig, hm)

fig
```
\end{examplefigure}

### With Mutating Plotting Functions

Mutating plotting functions work a bit differently with GridPositions.
First, it is checked if one - and only one - axis exists already at the given position.
If that's the case, that axis is plotted into.
If it's not the case, the function will error.

\begin{examplefigure}{}
```julia
using GLMakie
GLMakie.activate!() # hide
fig = Figure()

lines(fig[1, 1], 1.0..10, sin, color = :blue)
# this works because the previous command created an axis at fig[1, 1]
lines!(fig[1, 1], 1.0..10, cos, color = :red)

# the following line wouldn't work yet because no axis exists at fig[1, 2]
# lines!(fig[1, 2], 1.0..10, sin, color = :green)

fig[1, 2] = Axis(fig)
# now it works
lines!(fig[1, 2], 1.0..10, sin, color = :green)

# also works with nested grids
fig[2, 1:2][1, 3] = Axis(fig)
lines!(fig[2, 1:2][1, 3], 1.0..10, cos, color = :orange)

# but often it's more convenient to save an axis to reuse it
ax, _ = lines(fig[2, 1:2][1, 1:2], 1.0..10, sin, color = :black)
lines!(ax, 1.0..10, cos, color = :yellow)

fig
```
\end{examplefigure}
# Plot Recipes

Recipes allow you to extend `Makie` with your own custom types and plotting commands.

There are two types of recipes:

- _Type recipes_ define a simple mapping from a user defined type to an existing plot type
- _Full recipes_ define new custom plotting functions.

## Type recipes

Type recipes are mostly just conversions from one type or set of input argument types, yet unknown to Makie, to another which Makie can handle already.

This is the sequential logic by which conversions in Makie are attempted:

- Dispatch on `convert_arguments(::PlotType, args...)`
- If no matching method is found, determine a conversion trait via `conversion_trait(::PlotType)`
- Dispatch on `convert_arguments(::ConversionTrait, args...)`
- If no matching method is found, try to convert each single argument recursively with `convert_single_argument` until each type doesn't change anymore
- Dispatch on `convert_arguments(::PlotType, converted_args...)`
- Fail if no method was found

### Multiple Argument Conversion with `convert_arguments`

Plotting of a `Circle` for example can be defined via a conversion into a vector of points:

```julia
Makie.convert_arguments(x::Circle) = (decompose(Point2f, x),)
```

!!! warning
    `convert_arguments` must always return a Tuple.

You can restrict conversion to a subset of plot types, like only for scatter plots:

```julia
Makie.convert_arguments(P::Type{<:Scatter}, x::MyType) = convert_arguments(P, rand(10, 10))
```

Conversion traits make it easier to define behavior for a group of plot types that share the same trait. `PointBased` for example applies to `Scatter`, `Lines`, etc. Predefined are `NoConversion`, `PointBased`, `SurfaceLike` and `VolumeLike`.

```julia
Makie.convert_arguments(P::PointBased, x::MyType) = ...
```

Lastly, it is also possible to convert multiple arguments together.

```julia
Makie.convert_arguments(P::Type{<:Scatter}, x::MyType, y::MyOtherType) = ...
```

Optionally you may define the default plot type so that `plot(x::MyType)` will
use it directly:

```julia
plottype(::MyType) = Surface
```

### Single Argument Conversion with `convert_single_argument`

Some types which are unknown to Makie can be converted to other types, for which `convert_arguments` methods are available.
This is done with `convert_single_argument`.

For example, `AbstractArrays` with `Real`s and `missing`s can usually be safely converted to `Float32` arrays with `NaN`s instead of `missing`s.

The difference between `convert_single_argument` and `convert_arguments` with a single argument is that the former can be applied to any argument of any signature, while the latter only matches one-argument signatures.

## Full recipes with the `@recipe` macro

A full recipe comes in two parts. First is the plot type name, for example `MyPlot`, and then arguments and theme definition which are defined using the `@recipe` macro.

Second is at least one custom `plot!` method for `MyPlot` which creates an actual visualization using other existing plotting functions.

We use an example to show how this works:

```julia
@recipe(MyPlot, x, y, z) do scene
    Theme(
        plot_color = :red
    )
end
```

This macro expands to several things. Firstly a type definition:

```julia
const MyPlot{ArgTypes} = Combined{myplot, ArgTypes}
```

The type parameter of `Combined` contains the function `myplot` instead of e.g. a
symbol `MyPlot`. This way the mapping from `MyPlot` to `myplot` is safer and simpler.
The following signatures are automatically defined to make `MyPlot` nice to use:

```julia
myplot(args...; kw_args...) = ...
myplot!(args...; kw_args...) = ...
```

A specialization of `argument_names` is emitted if you have an argument list
`(x,y,z)` provided to the recipe macro:

```julia
argument_names(::Type{<: MyPlot}) = (:x, :y, :z)
```

This is optional but it will allow the use of `plot_object[:x]` to
fetch the first argument from the call
`plot_object = myplot(rand(10), rand(10), rand(10))`, for example.

Alternatively you can always fetch the `i`th argument using `plot_object[i]`,
and if you leave out the `(x,y,z)`, the default version of `argument_names`
will provide `plot_object[:arg1]` etc.

The theme given in the body of the `@recipe` invocation is inserted into a
specialization of `default_theme` which inserts the theme into any scene that
plots `MyPlot`:

```julia
function default_theme(scene, ::MyPlot)
    Theme(
        plot_color => :red
    )
end
```

As the second part of defining `MyPlot`, you should implement the actual
plotting of the `MyPlot` object by specializing `plot!`:

```julia
function plot!(myplot::MyPlot)
    # normal plotting code, building on any previously defined recipes
    # or atomic plotting operations, and adding to the combined `myplot`:
    lines!(myplot, rand(10), color = myplot[:plot_color])
    plot!(myplot, myplot[:x], myplot[:y])
    myplot
end
```

It's possible to add specializations here, depending on the argument _types_
supplied to `myplot`. For example, to specialize the behavior of `myplot(a)`
when `a` is a 3D array of floating point numbers:

```julia
const MyVolume = MyPlot{Tuple{<:AbstractArray{<: AbstractFloat, 3}}}
argument_names(::Type{<: MyVolume}) = (:volume,) # again, optional
function plot!(plot::MyVolume)
    # plot a volume with a colormap going from fully transparent to plot_color
    volume!(plot, plot[:volume], colormap = :transparent => plot[:plot_color])
    plot
end
```

## Example: Stock Chart

Let's say we want to visualize stock values with the classic open / close and low / high combinations.
In this example, we will create a special type to hold this information, and a recipe that can plot this type.

First, we make a struct to hold the stock's values for a given day:

```julia:eval-env
using CairoMakie
CairoMakie.activate!() # hide

struct StockValue{T<:Real}
    open::T
    close::T
    high::T
    low::T
end
```

Now we create a new plot type called `StockChart`.
The `do scene` closure is just a function that returns our default attributes, in this case they color stocks going down red, and stocks going up green.

```julia:eval-env
@recipe(StockChart) do scene
    Attributes(
        downcolor = :red,
        upcolor = :green,
    )
end
nothing # hide
```

Then we get to the meat of the recipe, which is actually creating a plot method.
We need to overload a specific method of `Makie.plot!` which as its argument has a subtype of our new `StockChart` plot type.
The type parameter of that type is a Tuple describing the argument types for which this method should work.

Note that the input arguments we receive inside the `plot!` method, which we can extract by indexing into the `StockChart`, are automatically converted to Observables by Makie.

This means that we must construct our plotting function in a dynamic way so that it will update itself whenever the input observables change.
This can be a bit trickier than recipes you might know from other plotting packages which produce mostly static plots.

```julia:eval-env
function Makie.plot!(
        sc::StockChart{<:Tuple{AbstractVector{<:Real}, AbstractVector{<:StockValue}}})

    # our first argument is an observable of parametric type AbstractVector{<:Real}
    times = sc[1]
    # our second argument is an observable of parametric type AbstractVector{<:StockValue}}
    stockvalues = sc[2]

    # we predefine a couple of observables for the linesegments
    # and barplots we need to draw
    # this is necessary because in Makie we want every recipe to be interactively updateable
    # and therefore need to connect the observable machinery to do so
    linesegs = Observable(Point2f[])
    bar_froms = Observable(Float32[])
    bar_tos = Observable(Float32[])
    colors = Observable(Bool[])

    # this helper function will update our observables
    # whenever `times` or `stockvalues` change
    function update_plot(times, stockvalues)
        colors[]

        # clear the vectors inside the observables
        empty!(linesegs[])
        empty!(bar_froms[])
        empty!(bar_tos[])
        empty!(colors[])

        # then refill them with our updated values
        for (t, s) in zip(times, stockvalues)
            push!(linesegs[], Point2f(t, s.low))
            push!(linesegs[], Point2f(t, s.high))
            push!(bar_froms[], s.open)
            push!(bar_tos[], s.close)
        end
        append!(colors[], [x.close > x.open for x in stockvalues])
        colors[] = colors[]
    end

    # connect `update_plot` so that it is called whenever `times`
    # or `stockvalues` change
    Makie.Observables.onany(update_plot, times, stockvalues)

    # then call it once manually with the first `times` and `stockvalues`
    # contents so we prepopulate all observables with correct values
    update_plot(times[], stockvalues[])

    # for the colors we just use a vector of booleans or 0s and 1s, which are
    # colored according to a 2-element colormap
    # we build this colormap out of our `downcolor` and `upcolor`
    # we give the observable element type `Any` so it will not error when we change
    # a color from a symbol like :red to a different type like RGBf(1, 0, 1)
    colormap = lift(Any, sc.downcolor, sc.upcolor) do dc, uc
        [dc, uc]
    end

    # in the last step we plot into our `sc` StockChart object, which means
    # that our new plot is just made out of two simpler recipes layered on
    # top of each other
    linesegments!(sc, linesegs, color = colors, colormap = colormap)
    barplot!(sc, times, bar_froms, fillto = bar_tos, color = colors, strokewidth = 0, colormap = colormap)

    # lastly we return the new StockChart
    sc
end
nothing # hide
```

Finally, let's try it out and plot some stocks:

\begin{examplefigure}{}
```julia
timestamps = 1:100

# we create some fake stock values in a way that looks pleasing later
startvalue = StockValue(0.0, 0.0, 0.0, 0.0)
stockvalues = foldl(timestamps[2:end], init = [startvalue]) do values, t
    open = last(values).close + 0.3 * randn()
    close = open + randn()
    high = max(open, close) + rand()
    low = min(open, close) - rand()
    push!(values, StockValue(
        open, close, high, low
    ))
end

# now we can use our new recipe
f = Figure()

stockchart(f[1, 1], timestamps, stockvalues)

# and let's try one where we change our default attributes
stockchart(f[2, 1], timestamps, stockvalues,
    downcolor = :purple, upcolor = :orange)
f
```
\end{examplefigure}

As a last example, lets pretend our stock data is coming in dynamically, and we want to create an animation out of it.
This is easy if we use observables as input arguments which we then update frame by frame:

```julia:stockchart_animation
timestamps = Observable(collect(1:100))
stocknode = Observable(stockvalues)

fig, ax, sc = stockchart(timestamps, stocknode)

record(fig, "stockchart_animation.mp4", 101:200,
        framerate = 30) do t
    # push a new timestamp without triggering the observable
    push!(timestamps[], t)

    # push a new StockValue without triggering the observable
    old = last(stocknode[])
    open = old.close + 0.3 * randn()
    close = open + randn()
    high = max(open, close) + rand()
    low = min(open, close) - rand()
    new = StockValue(open, close, high, low)
    push!(stocknode[], new)

    # now both timestamps and stocknode are synchronized
    # again and we can trigger one of them by assigning it to itself
    # to update the whole stockcharts plot for the new frame
    stocknode[] = stocknode[]
    # let's also update the axis limits because the plot will grow
    # to the right
    autolimits!(ax)
end
nothing # hide

using GLMakie # hide
GLMakie.activate!() # hide
```

\video{stockchart_animation}
# How layouts work

The goal of MakieLayout is that all elements placed in a scene fit into the
window, fill the available space, and are nicely aligned relative to each other.
This works by using `GridLayout` objects that determine how wide their rows and
columns should be given their content elements.

Content elements have inner widths and heights, as well as four protrusions, that tell
how far supporting content (like axis decorations) sticks out from the main part.
The protrusions are meant to stick into the gaps between grid cells, and not every
element has meaningful protrusions. They are mostly meant to allow for alignment
of axes along their spines.

Each element in a layout should have a couple of observables that support the layout
computations.
- Suggested bounding box
- Computed bounding box
- Auto-determined width and height
- Computed width and height
- Protrusions
- Size attributes
- Alignment attributes

### Suggested bounding box

This is the bounding box that is suggested to the element. Depending on the
settings of the element, it can choose to align perfectly with this bounding box
or, if its actual dimensions differ, how it should align inside that rectangle.
A small `Label` can for example be aligned top-left inside a big available suggested
bounding box.

### Computed bounding box

This is the bounding box of the element after it has received a suggested bounding
box and applied its own layout logic. This is the bounding box in which the elements
main area will be in the scene.

### Auto-determined width and height

Some elements can compute their own size, depending on their settings. `Label`,
for example, can compute the bounding box of its text. If an object has no specific
content, like an `Axis`, the auto-determined width or height will be `nothing`.

### Computed width and height

The computed width and height is the size that the element reports to a `GridLayout`
that it is a content element of. This can be different from the auto-size if the
object doesn't want its parent layout to know its auto-size. This is useful if
you don't want a column to shrink to the size of a `Label`, for example.

### Protrusions

These are four values that tell the `GridLayout` how much gap space is needed by
the element outside of the main element area. With an `Axis` that would be the
title at the top, y axis at the left side and x axis at the bottom in standard
configuration.

### Size attributes

The user can specify height and width of an element in different ways, which interact with the
suggested bounding box and the auto-determined size to compute the final size of the object and
also control how the layout responds to the element's size (used here for either width or height, respectively).

- `Fixed` or `Real`: The size is always fixed, no matter what the layout suggests. A `GridLayout` can auto-adjust column sizes to this size.
- `Relative`: The size is a fraction of the suggested size. A `GridLayout` can not auto-adjust column sizes to this size.
- `Auto`: The size is equal to the auto-determined size if it's not `nothing`. A `GridLayout` can auto-adjust to this size if it's not `nothing`.
- `nothing`: The size is equal to the suggested size. A `GridLayout` can not auto-adjust column sizes to this size.

For all sizes that a `GridLayout` can auto-adjust to, you can prohibit that by setting
`tellheight` or `tellwidth` of the element to `false`.

### Alignment attributes

The user can specify how an element should be aligned relative to its suggested
bounding box if it's not of the same size (in which case the alignment just has no effect on placement).
Currently, these values can be `:left`, `:right` or `:center` for horizontal alignment
and `:top`, `:bottom` and `:center` for vertical alignment.
# Figures

The `Figure` object contains a top-level `Scene` and a `GridLayout`, as well as a list of layoutables that have been placed into it, like `Axis`, `Colorbar`, `Slider`, `Legend`, etc.


## Creating a `Figure`

You can create a figure explicitly with the `Figure()` function, and set attributes of the underlying scene.
The most important one of which is the `resolution`.

```julia
f = Figure()
f = Figure(resolution = (600, 400))
```

A figure is also created implicitly when you use simple, non-mutating plotting commands like `plot()`, `scatter()`, `lines()`, etc.
Because these commands also create an axis for the plot to live in and the plot itself, they return a compound object `FigureAxisPlot`, which just stores these three parts.
To access the figure you can either destructure that object into its three parts or access the figure field directly.

```julia
figureaxisplot = scatter(rand(100, 2))
figure = figureaxisplot.figure

# destructuring syntax
figure, axis, plot = scatter(rand(100, 2))

# you can also ignore components
figure, = scatter(rand(100, 2))
```

You can pass arguments to the created figure in a dict-like object to the special `figure` keyword:

```julia
scatter(rand(100, 2), figure = (resolution = (600, 400),))
```

## Placing layoutables into a `Figure`

All layoutables take their parent figure as the first argument, then you can place them in the figure layout via indexing syntax.

```julia
f = Figure()
ax = f[1, 1] = Axis(f)
sl = f[2, 1] = Slider(f)
```

## GridPositions and GridSubpositions

The indexing syntax of `Figure` is implemented to work seamlessly with layouting.
If you index into the figure, a `GridPosition` object that stores this indexing operation is created.
This object can be used to plot a new axis into a certain layout position in the figure, for example like this:

\begin{examplefigure}{}
```julia
using CairoMakie
CairoMakie.activate!() # hide

f = Figure()
pos = f[1, 1]
scatter(pos, rand(100, 2))

pos2 = f[1, 2]
lines(pos2, cumsum(randn(100)))

# you don't have to store the position in a variable first, of course
heatmap(f[1, 3], randn(10, 10))

f
```
\end{examplefigure}


You can also index further into a `GridPosition`, which creates a `GridSubposition`.
With `GridSubposition`s you can describe positions in arbitrarily nested grid layouts.
Often, a desired plot layout can only be achieved with nesting, and repeatedly indexing makes this easy.

\begin{examplefigure}{}
```julia
using CairoMakie
CairoMakie.activate!() # hide

f = Figure()

f[1, 1] = Axis(f, title = "I'm not nested")
f[1, 2][1, 1] = Axis(f, title = "I'm nested")

# plotting into nested positions also works
heatmap(f[1, 2][2, 1], randn(20, 20))

f
```
\end{examplefigure}


All nested GridLayouts that don't exist yet, but are needed for a nested plotting call, are created in the background automatically.

!!! note
    The `GridLayout`s that are implicitly created when using `GridSubpositions` are not directly available in the return
    value for further manipulation. You can instead retrieve them after the fact with the `content` function, for example,
    as explained in the following section.


## Figure padding

You can change the amount of whitespace around the figure content with the keyword `figure_padding`.
This takes either a number for all four sides, or a tuple of four numbers for left, right, bottom, top.
You can also theme this setting with `set_theme!(figure_padding = 30)`, for example.

\begin{examplefigure}{}
```julia
using CairoMakie
CairoMakie.activate!() # hide

f = Figure(figure_padding = 1, backgroundcolor = :gray80)

Axis(f[1, 1])
scatter!(1:10)

f
```
\end{examplefigure}

## Retrieving Objects From A Figure

Sometimes users are surprised that indexing into a figure does not retrieve the object placed at that position.
This is because the `GridPosition` is needed for plotting, and returning content objects directly would take away that possibility.
Furthermore, a `GridLayout` can hold multiple objects at the same position, or have partially overlapping content,
so it's not well-defined what should be returned given a certain index.

To retrieve objects from a Figure you can instead use indexing plus the `contents` or `content` functions.
The `contents` function returns a Vector of all objects found at the given `GridPosition`.
You can use the `exact = true` keyword argument so that the position has to match exactly, otherwise objects
contained in that position are also returned.

```julia
f = Figure()
box = f[1:3, 1:2] = Box(f)
ax = f[1, 1] = Axis(f)

contents(f[1, 1]) == [ax]
contents(f[1:3, 1:2]) == [box, ax]
contents(f[1:3, 1:2], exact = true) == [box]
```

If you use `contents` on a `GridSubposition`, the `exact` keyword only refers to the lowest-level
grid layout, all upper levels have to match exactly.

```julia
f = Figure()
ax = f[1, 1][2, 3] = Axis(f)

contents(f[1, 1][2, 3]) == [ax]
contents(f[1:2, 1:2][2, 3]) == [] # the upper level has to match exactly
```

Often, you will expect only one object at a certain position and you want to work directly with it, without
retrieving it from the Vector returned by `contents`.
In that case, use the `content` function instead.
It works equivalently to `only(contents(pos, exact = true))`, so it errors if it can't return exactly one object
from an exact given position.

```julia
f = Figure()
ax = f[1, 1] = Axis(f)

contents(f[1, 1]) == [ax]
content(f[1, 1]) == ax
```
# Scenes

## What is a `Scene`?

`Scene`s are fundamental building blocks of Makie figures.
A Scene is like a container for `Plot`s and other `Scene`s.
`Scenes` have `Plot`s and `Subscenes` associated with them.
Every Scene also has a transformation, made up of _scale_, _translation_, and _rotation_.

!!! note
    Before the introduction of the `Figure` workflow, `Scene`s used to be the main container object which was returned from all plotting functions.
    Now, scenes are mostly an implementation detail for many users, unless they want to build custom solutions that go beyond what the default system offers.

A Scene's plots can be accessed via `scene.plots`.

A Scene's subscenes (also called children) can be accessed through `scene.children`.  This will return an Array of the `Scene`'s child scenes.  A child scene can be created by `childscene = Scene(parentscene)`.

Any `Scene` with an axis also has a `camera` associated with it; this can be accessed through `camera(scene)`, and its controls through `cameracontrols(scene)`.  More documentation about these is in the \myreflink{Cameras} section.

`Scene`'s also have configurable size/resolution. You can set the size in pixels by doing `Scene(resolution = (500, 500))`.

Any keyword argument given to the `Scene` will be propagated to its plots; therefore, you can set the palette or the colormap in the Scene itself.

## Subscenes

A subscene is no different than a normal Scene, except that it is linked to a "parent" Scene.  It inherits the transformations of the parent Scene, but can then be transformed independently of it.

## Scene Attributes

* `scene.clear = true`: Scenes are drawn parent first onto the same image. If `clear = true` for a (sub)scene it will clear the previously drawn things in its region to its `backgroundcolor`. Otherwise the plots in `scene` will be drawn on top and the backgroundcolor will be ignored. Note that this is not technically an attribute but just a field of `Scene`.
* `ssao = SSAO(bias = 0.025, blur=2, radius=0.5)`: Controls SSAO settings, see lighting documentation.
* `resolution = (800, 600)`: Sets the size of the created window if the scene is the root scene.

## Modifying A Scene

Makie offers mutation functions to scale, translate and rotate your Scenes on the fly.

{{doc translate!}}
{{doc rotate!}}
{{doc scale!}}

## Updating the Scene

When the Scene is changed, you may need to update several aspects of it.
Makie provides three main updating functions:

{{doc update_cam!}}

## Events

Scenes have several pre-created event "hooks" (through Observables) that you can handle.  These can be accessed through `scene.events`, which returns an \apilink{Events} struct.
# Frequently Asked Questions

## Installation Issues

We assume you are running Julia on the default system image without PackageCompiler.

### No `Scene` displayed or GLMakie fails to build

If `Makie` builds, but when plotting no window or plot is displayed, your backend may not have built correctly.
By default, Makie will try to use GLMakie as a backend, but if it does not build correctly for whatever reason, then scenes will not be displayed.
Ensure that your graphics card supports OpenGL; if it does not (old models, or relatively old integrated graphics cards), then you may want to consider CairoMakie.

## Plotting issues

### Dimensions too large

In general, plotting functions tend to plot whatever's given to them as a single texture.  This can lead to GL errors, or OpenGL failing silently.  To circumvent this, one can 'tile' the plots (i.e., assemble them piece-by-piece) to decrease the individual texture size.

#### 2d plots (heatmaps, images, etc.)

```julia
heatmap(rand(Float32, 24900, 26620))
```
may either fail with an error
```julia
   Error showing value of type Scene:
ERROR: glTexImage 2D: width too large. Width: 24900
[...]
```
or fail silently:

![untiled heatmap](https://user-images.githubusercontent.com/32143268/55675737-96357280-5894-11e9-9170-1ffd21f544cc.png)

Tiling the plot, as shown below, yields a correct image.

```julia
sc = Scene()
data = rand(Float32, 24900, 26620)
heatmap!(sc, 1:size(data, 1)÷2, 1:size(data, 2)÷2, data[1:end÷2, 1:end÷2])
heatmap!(sc, (size(data, 1)÷2 + 1):size(data, 1), 1:size(data, 2)÷2, data[(end÷2 + 1):end, 1:end÷2])
heatmap!(sc, 1:size(data, 1)÷2, (size(data, 2)÷2 + 1):size(data, 2), data[1:end÷2, (end÷2 + 1):end])
heatmap!(sc, (size(data, 1)÷2 + 1):size(data, 1), (size(data, 2)÷2 + 1):size(data, 2),
         data[(end÷2 + 1):end, (end÷2 + 1):end])
```
![tiled heatmap](https://user-images.githubusercontent.com/32143268/61105143-a3b35780-a496-11e9-83d1-bebe549aa593.png)

#### 3d plots (volumes)

The approach here is similar to that for the 2d plots, except that here there is a helpful function that gives the maximum texture size.
You can check the maximum texture size with:
```julia
using Makie, GLMakie, ModernGL
# simple plot to open a window (needs to be open for opengl)
display(scatter(rand(10)))
glGetIntegerv(GL_MAX_3D_TEXTURE_SIZE)
```
and then just split the volume:
```julia
vol = rand(506, 720, 1440)
ranges = (1:256, 1:256, 1:256)
scene = volume(ranges..., vol[ranges...])
for i in 1:3
    global ranges
    ranges = ntuple(3) do j
        s = j == i ? last(ranges[j]) : 1
        e = j == i ? size(vol, j) : last(ranges[j])
        s:e
    end
    volume!(ranges..., vol[ranges...])
end
scene
```

## General issues

### My font doesn't work!

If `Makie` can't find your font, you can do two things:

1) Check that the name matches and that the font is in one of the directories in:

    - `using FreeTypeAbstraction; FreeTypeAbstraction.valid_fontpaths`

2) You can add a custom font path via the environment variable:

    - `ENV["FREETYPE_ABSTRACTION_FONT_PATH"] = "/path/to/your/fonts"`

3) Specify the path to the font; instead of `font = "Noto"`, you could write `joindir(homedir(), "Noto.ttf")` or something.


## Layout Issues

### Elements are squashed into the lower left corner

Layoutable elements require a bounding box that they align themselves to. If you
place such an element in a layout, the bounding box is controlled by that layout.
If you forget to put an element in a layout, it will have its default bounding box
of `BBox(0, 100, 0, 100)` which ends up being in the lower left corner. You can
also choose to specify a bounding box manually if you need more control.

\begin{examplefigure}{svg = true}
```julia
using CairoMakie
CairoMakie.activate!() # hide

f = Figure()

ax1 = Axis(f, title = "Squashed")
ax2 = Axis(f[1, 1], title = "Placed in Layout")
ax3 = Axis(f, bbox = BBox(200, 600, 100, 500),
  title = "Placed at BBox(200, 600, 100, 500)")

f
```
\end{examplefigure}


### Columns or rows are shrunk to the size of Text or another element

Columns or rows that have size `Auto(true)` try to determine the width or height of all
single-spanned elements that are placed in them, and if any elements report their
size the row or column will shrink to the maximum reported size. This is so smaller
elements with a known size take as little space as needed. But if there is other
content in the row that should take more space, you can give the offending element
the attribute `tellheight = false` or `tellwidth = false`. This way, its own size
can be determined automatically, but
it doesn't report it to the row or column of the layout. Alternatively, you can set the size
of that row or column to `Auto(false)` (or any other value than `Auto(true)`).

\begin{examplefigure}{svg = true}
```julia
using CairoMakie

f = Figure()

Axis(f[1, 1], title = "Shrunk")
Axis(f[2, 1], title = "Expanded")
Label(f[1, 2], "tellheight = true", tellheight = true)
Label(f[2, 2], "tellheight = false", tellheight = false)

f
```
\end{examplefigure}


### The Figure content does not fit the Figure

`GridLayout`s work by fitting all their child content into the space that is available to them.
Therefore, the `Figure` size determines how the layout is solved, but the layout does not influence the `Figure` size.

This works well when all content is adjustable in width and height, such as an `Axis` that can shrink or grow as needed.
But it is also possible to constrain elements or rows/columns in width, height or aspect.
And if too many elements have such constraints, it's not possible any longer to fit them all into the given `Figure` size, without leaving whitespace or clipping them at the borders.

If this is the case, you can use the function `resize_to_layout!`, which determines the actual size of the main `GridLayout` given its content, and resizes the `Figure` to fit.

Here is an example, where all `Axis` objects are given fixed widths and heights.
There are not enough degrees of freedom for the layout algorithm to fit everything nicely into the `Figure`:

\begin{examplefigure}{svg = true}
```julia
using CairoMakie

set_theme!(backgroundcolor = :gray90)

f = Figure(resolution = (800, 600))

for i in 1:3, j in 1:3
    ax = Axis(f[i, j], title = "$i, $j", width = 100, height = 100)
    i < 3 && hidexdecorations!(ax, grid = false)
    j > 1 && hideydecorations!(ax, grid = false)
end

Colorbar(f[1:3, 4])

f
```
\end{examplefigure}


As you can see, there's empty space on all four sides, because there are no flexible objects that could fill it.

But once we run `resize_to_layout!`, the `Figure` assumes the appropriate size for our axes:

\begin{examplefigure}{svg = true}
```julia
resize_to_layout!(f)
set_theme!() # hide
f
```
\end{examplefigure}# Observables & Interaction

Interaction and animations in Makie are handled using [`Observables.jl`](https://juliagizmos.github.io/Observables.jl/stable/).
An `Observable` is a container object whose stored value you can update interactively.
You can create functions that are executed whenever an observable changes.
You can also create observables whose values are updated whenever other observables change.
This way you can easily build dynamic and interactive visualizations.

On this page you will learn how the `Observable`s pipeline and the event-based interaction system work. Besides this, there is also a video tutorial on how to make interactive visualizations (or animations) with Makie.jl and the `Observable` system:

~~~<a class="boxlink" href="https://www.youtube.com/watch?v=L-gyDvhjzGQ">~~~
@@title Animations & Interaction @@
@@box-content
    @@description
    How to create animations and interactive applications in Makie.
    @@
    ~~~
    <div class="youtube-container">
    <iframe width="560" height="315" src="https://www.youtube.com/embed/L-gyDvhjzGQ?controls=0" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe>
    </div>
    ~~~
@@
~~~</a>~~~

## The `Observable` structure

A `Observable` is an object that allows its value to be updated interactively.
Let's start by creating one:

```julia:code1
using GLMakie, Makie

x = Observable(0.0)
```
\show{code1}

Each `Observable` has a type parameter, which determines what kind of objects it can store.
If you create one like we did above, the type parameter will be the type of the argument.
Keep in mind that sometimes you want a wider parametric type because you intend to update the `Observable` later with objects of different types.
You could for example write:

```julia
x2 = Observable{Real}(0.0)
x3 = Observable{Any}(0.0)
```

This is often the case when dealing with attributes that can come in different forms.
For example, a color could be `:red` or `RGB(1,0,0)`.

## Triggering A Change

You change the value of a Observable with empty index notation:

```julia:code2
x[] = 3.34
nothing # hide
```
\show{code2}

This was not particularly interesting.
But Observables allow you to register functions that are executed whenever the Observable's content is changed.

One such function is `on`. Let's register something on our Observable `x` and change `x`'s value:

```julia:code3
on(x) do x
    println("New value of x is $x")
end

x[] = 5.0
nothing # hide
```
\show{code3}


!!! note
    All registered functions in a `Observable` are executed synchronously in the order of registration.
    This means that if you change two Observables after one another, all effects of the first change will happen before the second change.

There are two ways to access the value of a `Observable`.
You can use the indexing syntax or the `to_value` function:

```julia
value = x[]
value = to_value(x)
```

The advantage of using `to_value` is that you can use it in situations where you could either be dealing with Observables or normal values. In the latter case, `to_value` just returns the original value, like `identity`.

## Chaining `Observable`s With `lift`

You can create a Observable depending on another Observable using \apilink{lift}.
The first argument of `lift` must be a function that computes the value of the output Observable given the values of the input Observables.

```julia:code4
f(x) = x^2
y = lift(f, x)
```
\show{code4}

Now, whenever `x` changes, the derived `Observable` `y` will immediately hold the value `f(x)`.
In turn, `y`'s change could trigger the update of other observables, if any have been connected.
Let's connect one more observable and update x:

```julia:code5
z = lift(y) do y
    -y
end

x[] = 10.0

@show x[]
@show y[]
@show z[]
nothing # hide
```
\show{code5}

If `x` changes, so does `y` and then `z`.

Note, though, that changing `y` does not change `x`.
There is no guarantee that chained Observables are always synchronized, because they
can be mutated in different places, even sidestepping the change trigger mechanism.

```julia:code6
y[] = 20.0

@show x[]
@show y[]
@show z[]
nothing # hide
```
\show{code6}


## Shorthand Macro For `lift`

When using \apilink{lift}, it can be tedious to reference each participating `Observable`
at least three times, once as an argument to `lift`, once as an argument to the closure that
is the first argument, and at least once inside the closure:

```julia
x = Observable(rand(100))
y = Observable(rand(100))
z = lift((x, y) -> x .+ y, x, y)
```

To circumvent this, you can use the `@lift` macro. You simply write the operation
you want to do with the lifted `Observable`s and prepend each `Observable` variable
with a dollar sign \$. The macro will lift every Observable variable it finds and wrap
the whole expression in a closure. The equivalent to the above statement using `@lift` is:

```julia
z = @lift($x .+ $y)
```

This also works with multiline statements and tuple or array indexing:

```julia
multiline_node = @lift begin
    a = $x[1:50] .* $y[51:100]
    b = sum($z)
    a .- b
end
```

If the Observable you want to reference is the result of some expression, just use `$` with parentheses around that expression.

```julia
container = (x = Observable(1), y = Observable(2))

@lift($(container.x) + $(container.y))
```

## Problems With Synchronous Updates

One very common problem with a pipeline based on multiple observables is that you can only change observables one by one.
Theoretically, each observable change triggers its listeners immediately.
If a function depends on two or more observables, changing one right after the other would trigger it multiple times, which is often not what you want.

Here's an example where we define two Observables and lift a third one from them:

```julia
xs = Observable(1:10)
ys = Observable(rand(10))

zs = @lift($xs .+ $ys)
```

Now let's update both `xs` and `ys`:

```julia
xs[] = 2:11
ys[] = rand(10)
```

We just triggered `zs` twice, even though we really only intended one data update.
But this double triggering is only part of the problem.

Both `xs` and `ys` in this example had length 10, so they could still be added without a problem.
If we want to append values to xs and ys, the moment we change the length of one of them, the function underlying `zs` will error because of a shape mismatch.
Sometimes the only way to fix this situation, is to mutate the content of one observable without triggering its listeners, then triggering the second one.

```julia
xs.val = 1:11 # mutate without triggering listeners
ys[] = rand(11) # trigger listeners of ys (in this case the same as xs)
```

Use this technique sparingly, as it increases the complexity of your code and can make reasoning about it more difficult.
It also only works if you can still trigger all listeners correctly.
For example, if another observable listened only to `xs`, we wouldn't have updated it correctly in the above workaround.
Often, you can avoid length change problems by using arrays of containers like `Point2f` or `Vec3f` instead of synchronizing two or three observables of single element vectors manually.
# Basic transparency

To make a plot transparent you need to add an alpha value to its `color` or `colormap`.

\begin{examplefigure}{}
```julia
using CairoMakie, FileIO
CairoMakie.activate!() # hide

# color
fig, ax, p = image(0..11, -1..11, rotr90(FileIO.load(Makie.assetpath("cow.png"))))
scatter!(ax, 1:10,fill(10, 10), markersize = 30, color = :red)
scatter!(ax, 1:10, fill(9, 10), markersize = 30, color = (:red, 0.5))
scatter!(ax, 1:10, fill(8, 10), markersize = 30, color = RGBf(0.8, 0.6, 0.1))
scatter!(ax, 1:10, fill(7, 10), markersize = 30, color = RGBAf(0.8, 0.6, 0.1, 0.5))

# colormap
scatter!(ax, 1:10, fill(5, 10), markersize = 30, color = 1:10, colormap = :viridis)
scatter!(ax, 1:10, fill(4, 10), markersize = 30, color = 1:10, colormap = (:viridis, 0.5),)
scatter!(ax, 1:10, fill(3, 10), markersize = 30, color = 1:10, colormap = (:red, :orange),)
scatter!(ax, 1:10, fill(2, 10), markersize = 30, color = 1:10, colormap = ((:red, 0.5), (:orange, 0.5)))
cm = [RGBf(x^2, 1 - x^2, 0.2) for x in range(0, 1, length=100)]
scatter!(ax, 1:10, fill(1, 10), markersize = 30, color = 1:10, colormap = cm)
cm = [RGBAf(x^2, 1 - x^2, 0.2, 0.5) for x in range(0, 1, length=100)]
scatter!(ax, 1:10, fill(0, 10), markersize = 30, color = 1:10, colormap = cm)
fig
```
\end{examplefigure}


# Details and Problems with transparency

The color generated from two overlapping transparent objects depends on their order. Consider for example a red and blue marker with the same level of transparency. If the blue marker is in front we expect a more blue color where they overlap. If the red one is in front we expect a more red color. 

\begin{examplefigure}{}
```julia
using CairoMakie
CairoMakie.activate!() # hide

scene = Scene(resolution = (400, 275))
campixel!(scene)
scatter!(
    scene, [100, 200, 300], [100, 100, 100], 
    color = [RGBAf(1,0,0,0.5), RGBAf(0,0,1,0.5), RGBAf(1,0,0,0.5)], 
    markersize=150
)
scatter!(scene, Point2f(150, 175), color = (:green, 0.5), markersize=150)
p = scatter!(scene, Point2f(250, 175), color = (:green, 0.5), markersize=150)
translate!(p, 0, 0, -1)
scene
```
\end{examplefigure}

The graphic above follows three rules in terms of transparency:

1. If two plots are at different z-levels, the one with the higher level will be in front of the lower z-level. (The  green circle on the right is behind all other plots.)
2. The drawing order of plots at the same z-level matches their creation order. (The red and blue circles are behind the left green circle.)
3. Plot elements are drawn in order. (The left red circle is behind the middle blue circle which is behind the right red circle.)

The first rule follows from explicit sorting of plots. It is only done in 2D because a plot can have variable depth in 3D. The second and third rules apply in both cases. They will however frequently generate the wrong results in 3D. Take for example two planes rotated to have a varying depth value:

\begin{examplefigure}{}
```julia
using CairoMakie
CairoMakie.activate!() # hide

fig = Figure()
ax = LScene(fig[1, 1], show_axis=false)
p1 = mesh!(ax, Rect2f(-1.5, -1, 3, 3), color = (:red, 0.5), shading = false)
p2 = mesh!(ax, Rect2f(-1.5, -2, 3, 3), color = (:blue, 0.5), shading = false)
rotate!(p1, Vec3f(0, 1, 0), 0.1)
rotate!(p2, Vec3f(0, 1, 0), -0.1)
fig
```
\end{examplefigure}

\begin{examplefigure}{}
```julia
using GLMakie
GLMakie.activate!() # hide

fig = Figure()
ax = LScene(fig[1, 1], show_axis=false)
p1 = mesh!(ax, Rect2f(-1.5, -1, 3, 3), color = (:red, 0.5), shading = false)
p2 = mesh!(ax, Rect2f(-1.5, -2, 3, 3), color = (:blue, 0.5), shading = false)
rotate!(p1, Vec3f(0, 1, 0), 0.1)
rotate!(p2, Vec3f(0, 1, 0), -0.1)
fig
```
\end{examplefigure}

Both backends handle this wrong. CairoMakie seems to ignore depth and just draws the planes in plotting order. This isn't quite true - CairoMakie does consider depth on a per-plot and in some cases on a per-element basis (e.g. triangles in a 3D mesh). But it can't handle depth on a per pixel level. 

GLMakie on the other hand can handle depth on a per-pixel level, as evident by the correct order shown above. The problem with transparency here is that the order of colors applied to a pixel is not known a priori. GLMakie will draw the red plane first and record depth values for each pixel. Then it will draw the blue plane if it's in front of the other. Solving this exactly would require collecting colors and depth values per pixel, sorting them and then blending them in order. This would be very expensive and is therefore rarely done. 


## Order independent transparency

GLMakie implements an approximate scheme for blending transparent colors - [Order Independent Transparency](https://jcgt.org/published/0002/02/09/) (OIT). Instead of using the usual order dependent blending `alpha * color + (1 - alpha) * background_color` it uses a weighted sum with weights based on depth and alpha. You can turn on OIT by setting `transparency = true` for a given plot.

\begin{examplefigure}{}
```julia
using GLMakie
GLMakie.activate!() # hide

fig = Figure()
ax = LScene(fig[1, 1], show_axis=false)
p1 = mesh!(ax, Rect2f(-2, -2, 4, 4), color = (:red, 0.5), shading = false, transparency = true)
p2 = mesh!(ax, Rect2f(-2, -2, 4, 4), color = (:blue, 0.5), shading = false, transparency = true)
p3 = mesh!(ax, Rect2f(-2, -2, 4, 4), color = (:red, 0.5), shading = false, transparency = true)
for (dz, p) in zip((-1, 0, 1), (p1, p2, p3))
    translate!(p, 0, 0, dz)
end
fig
```
\end{examplefigure}

Being an approximate scheme OIT has some strengths and weaknesses. There are two significant drawbacks of OIT:
1. Blending always happens - even if a fully opaque color (alpha = 1) should hide another.
2. Blending isn't sharp - when two colors with the same alpha value are blended at similar depth values their output color will be similar.

\begin{examplefigure}{}
```julia
using GLMakie
GLMakie.activate!() # hide

fig = Figure(resolution = (800, 400))
ax1 = LScene(fig[1, 1], show_axis=false)
p1 = mesh!(ax1, Rect2f(-2, -2, 4, 4), color = :red, shading = false, transparency = true)
p2 = mesh!(ax1, Rect2f(-2, -2, 4, 4), color = :blue, shading = false, transparency = true)
p3 = mesh!(ax1, Rect2f(-2, -2, 4, 4), color = :red, shading = false, transparency = true)
for (dz, p) in zip((-1, 0, 1), (p1, p2, p3))
    translate!(p, 0, 0, dz)
end

ax2 = LScene(fig[1, 2], show_axis=false)
p1 = mesh!(ax2, Rect2f(-1.5, -1, 3, 3), color = (:red, 0.5), shading = false, transparency=true)
p2 = mesh!(ax2, Rect2f(-1.5, -2, 3, 3), color = (:blue, 0.5), shading = false, transparency=true)
rotate!(p1, Vec3f(0, 1, 0), 0.1)
rotate!(p2, Vec3f(0, 1, 0), -0.1)
fig
```
\end{examplefigure}

Note that you can mix opaque `transparency = false` plots with transparent OIT plots without problems. So the first issue is not really an issue for truly opaque plots but rather close to opaque plots.

Another problem you may run into is that part of your plot becomes black or white. This is a result of floats becoming infinite during the OIT calculation because of a large number of transparent colors being added close to the camera. You can fix this in two way:

1. Move the `near` clipping plane closer to the camera. For an `ax::LScene` this can be done by adjusting `ax.scene.camera_controls.near[]` closer to 0 and calling `update_cam!(ax.scene, ax.scene.camera_controls)`.
2. Reducing `GLMakie.transparency_weight_scale[]` which controls the maximum weight given to a transparent color. The default is `1000f0`.

Adjusting the latter may also help with sharpness.@def order = 99

# API reference

{{ api_reference }}
# Exporting a Figure with physical dimensions

Makie currently uses a unitless approach for specifying Figure resolution, font sizes, line widths, etc. The dimensions of the final result depend on which export format you use, and what settings you specify when saving. Unitless means that `resolution = (800, 600)` doesn't mean pixels, or mm, or cm. Before saving or displaying, these are just numbers.

GLMakie and WGLMakie can only export bitmaps (png files). CairoMakie can export both bitmaps and vector graphics (svg and pdf). These two file types have fundamentally different concepts of size.

## Bitmaps

Bitmaps always have a fixed resolution in pixels. Pixels are not a physical dimension. Monitor pixels have a physical size, and printer dots have one as well, but not the actual image. An image only has a physical size if you decide for a mapping from pixels to physical dimensions. That is usually the `dpi` value (dots per inch).

This value comes from printing and tells you how many printer dots fit into one inch. If we say that one printer dot corresponds to one pixel, then we can use `dpi` to convert from pixels to inch. (Pixels are not actually dots, but this interpretation is widely used. The actual metric to go from pixels to inches is `ppi` or pixels per inch.)

So if you need a plot at 4 x 3 inches and 400 dpi, you could set the figure size like this:

```!
size_in_inches = (4, 3)
dpi = 400
size_in_pixels = size_in_inches .* dpi
```

This means, if you export a bitmap, it doesn't have a physical size per se, but you can choose one later by deciding at what dpi you want to print the image. If you place a bitmap into a LaTeX document, for example, you can choose the size yourself, but just have to take care that the resulting dpi are high enough for your purposes, and that the document doesn't look blurry when printed.

When you save a `Figure` as a bitmap with GLMakie or WGLMakie, the unitless resolution can be interpreted as pixel size. If you save a bitmap with CairoMakie, you additionally have the option to use a scaling factor that decides the mapping from unitless dimensions to pixels. This is done with the `px_per_unit` keyword argument.

```julia
f = Figure(resolution = (800, 600))
# in GLMakie or WGLMakie
save("figure.png", f) # output size = 800 x 600 pixels
# in CairoMakie
save("figure.png", f) # output size = 800 x 600 pixels
save("figure.png", f, px_per_unit = 2) # output size = 1600 x 1200 pixels
```

## Vector graphics

Vector graphics don't have a resolution like bitmap images. They are collections of mathematical descriptions of lines and curves, and these can be arbitrarily scaled up and down. Vector graphics do have physical dimensions in that their content size is usually specified in `pt` which has a direct mapping to inch by convention, as 1 inch is equivalent to 72 pt.

When you export vector graphics with CairoMakie, you can control the mapping from unitless size to size in pt with the `pt_per_unit` keyword argument. This is by default set to `0.75`. The reason for this is that in web contexts, by convention 1 px is equivalent to 0.75 pt. So if you use Pluto or Jupyter notebooks and display a figure once as a bitmap and once as a vector graphic, they will have the same apparent size.

If you want to save a figure for a publication, you usually want to fit font sizes to the rest of the document, and adjust the size to what the journal expects.

Font sizes are usually given in pt, figure sizes in inches.

So if you need a 4 x 3 inches figure and your font size is 12 pt, you should set up and save your figure like this:

```julia
size_inches = (4, 3)
size_pt = 72 .* size_inches
f = Figure(resolution = size_pt, fontsize = 12)
save("figure.pdf", f, pt_per_unit = 1)
```# Cameras

A `Camera` is simply a viewport through which the Scene is visualized.  `Makie` offers 2D and 3D projections, and 2D plots can be projected in 3D!

To specify the camera you want to use for your Scene, you can set the `camera` attribute.  Currently, we offer four types of camera:

\apilink{campixel!}
\apilink{cam2d!}
`cam3d!`
`cam3d_cad!`

which will mutate the camera of the Scene into the specified type.

## Pixel Camera

The pixel camera (\apilink{campixel!(scene)}) projects the scene in pixel space, i.e. each integer step in the displayed data will correspond to one pixel. There are no controls for this camera. The clipping limits are set to `(-10_000, 10_000)`.

## 2D Camera

The 2D camera (\apilink{cam2d!(scene)}) uses an orthographic projection with a fixed rotation and aspect ratio. You can set the following attributes via keyword arguments in `cam2d!` or by accessing the camera struct `cam = cameracontrols(scene)`:

- `zoomspeed = 0.10f0` sets the speed of mouse wheel zooms.
- `zoombutton = nothing` sets an additional key that needs to be pressed in order to zoom. Defaults to no key.
- `panbutton = Mouse.right` sets the mouse button that needs to be pressed to translate the view.
- `selectionbutton = (Keyboard.space, Mouse.left)` sets a set of buttons that need to be pressed to perform rectangle zooms.

Note that this camera is not used by MakieLayout `Axis`. It is used, by default, for 2D `LScene`s and `Scene`s.

## 3D Camera

{{doc Camera3D}}

## General Remarks

To force a plot to be visualized in 3D, you can set the limits to have a nonzero \(z\)-axis interval, or ensure that a 3D camera type is used.
For example, you could pass the keyword argument `limits = Rect([0,0,0],[1,1,1])`, or `camera = cam3d!`.

Often, when modifying the Scene, the camera can get "out of sync" with the Scene. To fix this, you can call the \apilink{update_cam!} function on the Scene.

Buttons passed to the 2D and 3D camera are forwarded to `ispressed`. As such you can pass `false` to disable an interaction, `true` to ignore a modifier, any button, collection of buttons or even logical expressions of buttons. See the events documentation for more details.
# Lighting

For 3D scenes, `GLMakie` offers several attributes to control the lighting of the material.

- `ambient::Vec3f`: Objects should never be completely dark; we use an ambient light to simulate background lighting, and give the object some color. Each element of the vector represents the intensity of color in R, G or B respectively.
- `diffuse::Vec3f`: Simulates the directional impact which the light source has on the plot object. This is the most visually significant component of the lighting model; the more a part of an object faces the light source, the brighter it becomes. Each element of the vector represents the intensity of color in R, G or B respectively.
- `specular::Vec3f`: Simulates the bright spot of a light that appears on shiny objects. Specular highlights are more inclined to the color of the light than the color of the object. Each element of the vector represents the intensity of color in R, G or B respectively.
- `shininess::Float32`: Controls the shininess of the object. Higher shininess reduces the size of the highlight, and makes it sharper. This value must be positive.
- `lightposition::Vec3f`: The location of the main light source; by default, the light source is at the location of the camera.

You can find more information on how these were implemented [here](https://learnopengl.com/Lighting/Basic-Lighting).
Some usage examples can be found in the [RPRMakie examples](https://makie.juliaplots.org/stable/documentation/backends/rprmakie/) and in the [examples](https://makie.juliaplots.org/stable/documentation/lighting/#examples).

## SSAO

GLMakie also implements [_screen-space ambient occlusion_](https://learnopengl.com/Advanced-Lighting/SSAO), which is an algorithm to more accurately simulate the scattering of light. There are a couple of controllable scene attributes nested within the `SSAO` toplevel attribute:

- `radius` sets the range of SSAO. You may want to scale this up or
  down depending on the limits of your coordinate system
- `bias` sets the minimum difference in depth required for a pixel to
  be occluded. Increasing this will typically make the occlusion
  effect stronger.
- `blur` sets the (pixel) range of the blur applied to the occlusion texture.
  The texture contains a (random) pattern, which is washed out by
  blurring. Small `blur` will be faster, sharper and more patterned.
  Large `blur` will be slower and smoother. Typically `blur = 2` is
  a good compromise.

!!! note
    The SSAO postprocessor is turned off by default to save on resources. To turn it on, set `GLMakie.enable_SSAO[] = true`, close any existing GLMakie window and reopen it.

## Matcap

A matcap (material capture) is a texture which is applied based on the normals of a given mesh. They typically include complex materials and lighting and offer a cheap way to apply those to any mesh. You may pass a matcap via the `matcap` attribute of a `mesh`, `meshscatter` or `surface` plot. Setting `shading = false` is suggested. You can find a lot matcaps [here](https://github.com/nidorx/matcaps).

## Examples

\begin{showhtml}{}
```julia
using JSServe
Page(exportable=true, offline=true)
```
\end{showhtml}

\begin{showhtml}{}
```julia
using WGLMakie
using JSServe

WGLMakie.activate!() # hide
xs = -10:0.1:10
ys = -10:0.1:10
zs = [10 * (cos(x) * cos(y)) * (.1 + exp(-(x^2 + y^2 + 1)/10)) for x in xs, y in ys]

fig, ax, pl = surface(xs, ys, zs, colormap = (:white, :white),

    # Light comes from (0, 0, 15), i.e the sphere
    axis = (
        # Light comes from (0, 0, 15), i.e the sphere
        lightposition = Vec3f(0, 0, 15),
        # base light of the plot only illuminates red colors
        ambient = RGBf(0.3, 0, 0),
    ),
    # light from source (sphere) illuminates yellow colors
    diffuse = Vec3f(0.4, 0.4, 0),
    # reflections illuminate blue colors
    specular = Vec3f(0, 0, 1.0),
    # Reflections are sharp
    shininess = 128f0,
    figure = (resolution=(1000, 800),)
)
mesh!(ax, Sphere(Point3f(0, 0, 15), 1f0), color=RGBf(1, 0.7, 0.3))

app = JSServe.App() do session
    light_rotation = JSServe.Slider(1:360)
    shininess = JSServe.Slider(1:128)

    pointlight = ax.scene.lights[1]
    ambient = ax.scene.lights[2]
    on(shininess) do value
        pl.shininess = value
    end
    on(light_rotation) do degree
        r = deg2rad(degree)
        pointlight.position[] = Vec3f(sin(r)*10, cos(r)*10, 15)
    end
    JSServe.record_states(session, DOM.div(light_rotation, shininess, fig))
end
app
```
\end{showhtml}

\begin{examplefigure}{}
```julia
using GLMakie
GLMakie.activate!() # hide
GLMakie.enable_SSAO[] = true
close(GLMakie.global_gl_screen()) # close any open screen

fig = Figure()
ssao = Makie.SSAO(radius = 5.0, blur = 3)
ax = LScene(fig[1, 1], scenekw = (ssao=ssao,))
# SSAO attributes are per scene
ax.scene.ssao.bias[] = 0.025

box = Rect3(Point3f(-0.5), Vec3f(1))
positions = [Point3f(x, y, rand()) for x in -5:5 for y in -5:5]
meshscatter!(ax, positions, marker=box, markersize=1, color=:lightblue, ssao=true)
fig
```
\end{examplefigure}

```julia:disable-ssao
GLMakie.enable_SSAO[] = false # hide
close(GLMakie.global_gl_screen()) # hide
```

\begin{examplefigure}{}
```julia
using FileIO
using GLMakie
GLMakie.activate!() # hide
catmesh = FileIO.load(assetpath("cat.obj"))
gold = FileIO.load(download("https://raw.githubusercontent.com/nidorx/matcaps/master/1024/E6BF3C_5A4719_977726_FCFC82.png"))

mesh(catmesh, matcap=gold, shading=false)
```
\end{examplefigure}
# Colors

## Named colors
Named colors in Makie.jl (e.g., `:blue`) are parsed using [Colors.jl](https://juliagraphics.github.io/Colors.jl/stable/constructionandconversion/#Color-Parsing) and thus have a large array of possibilities under CSS specifications. You can find a plotted table of all possible names [in this page](https://juliagraphics.github.io/Colors.jl/stable/namedcolors/).

## Colormaps
{{colorschemes}}
# Predefined themes

Makie has a few predefined themes.
Here you can see the same example figure with these different themes applied.

~~~
<input id="hidecode" class="hidecode" type="checkbox">
~~~

```julia:demofigure
using CairoMakie
using Random

CairoMakie.activate!(type = "svg") # hide

function demofigure()
    Random.seed!(2)

    f = Figure()
    ax = Axis(f[1, 1],
        title = "measurements",
        xlabel = "time (s)",
        ylabel = "amplitude")

    labels = ["alpha", "beta", "gamma", "delta", "epsilon", "zeta"]
    for i in 1:6
        y = cumsum(randn(10)) .* (isodd(i) ? 1 : -1)
        lines!(y, label = labels[i])
        scatter!(y, label = labels[i])
    end

    Legend(f[1, 2], ax, "legend", merge = true)

    Axis3(f[1, 3],
        viewmode = :stretch,
        zlabeloffset = 40,
        title = "sinusoid")

    s = surface!(0:0.5:10, 0:0.5:10, (x, y) -> sqrt(x * y) + sin(1.5x))

    Colorbar(f[1, 4], s, label = "intensity")

    ax = Axis(f[2, 1:2],
        title = "different species",
        xlabel = "height (m)",
        ylabel = "density",)
    for i in 1:6
        y = randn(200) .+ 2i
        density!(y)
    end
    tightlimits!(ax, Bottom())
    xlims!(ax, -1, 15)

    Axis(f[2, 3:4],
        title = "stock performance",
        xticks = (1:6, labels),
        xlabel = "company",
        ylabel = "gain (\$)",
        xticklabelrotation = pi/6)
    for i in 1:6
        data = randn(1)
        barplot!([i], data)
        rangebars!([i], data .- 0.2, data .+ 0.2)
    end

    f
end
```

~~~
<label for="hidecode" class="hidecode"></label>
~~~

## Default theme

\begin{examplefigure}{}
```julia
demofigure()
```
\end{examplefigure}

## theme_ggplot2

\begin{examplefigure}{}
```julia
with_theme(demofigure, theme_ggplot2())
```
\end{examplefigure}

## theme_minimal

\begin{examplefigure}{}
```julia
with_theme(demofigure, theme_minimal())
```
\end{examplefigure}

## theme_black

\begin{examplefigure}{}
```julia
with_theme(demofigure, theme_black())
```
\end{examplefigure}

## theme_light

\begin{examplefigure}{}
```julia
with_theme(demofigure, theme_light())
```
\end{examplefigure}

## theme_dark

\begin{examplefigure}{}
```julia
with_theme(demofigure, theme_dark())
```
\end{examplefigure}
# CairoMakie

[CairoMakie](https://github.com/JuliaPlots/Makie.jl/tree/master/CairoMakie) uses Cairo.jl to draw vector graphics to SVG and PDF.
You should use it if you want to achieve the highest-quality plots for publications, as the rendering process of the GL backends works via bitmaps and is geared more towards speed than pixel-perfection.

### Special CairoMakie Properties

#### Inline Plot Type

You can choose the type of plot that is displayed inline in, e.g., VSCode, Pluto.jl, or any other environment, by setting it via the `activate!` function.

```julia
CairoMakie.activate!(type = "png")
CairoMakie.activate!(type = "svg")
```

#### Resolution Scaling

When you save a CairoMakie figure, you can change the mapping from figure resolution to pixels (when saving to png) or points (when saving to svg or pdf).
This way you can easily scale the resulting image up or down without having to change any plot element sizes.

Just specify `pt_per_unit` when saving vector formats and `px_per_unit` when saving pngs.
`px_per_unit` defaults to 1 and `pt_per_unit` defaults to 0.75.
When embedding svgs in websites, `1px` is equivalent to `0.75pt`.
This means that by default, saving a png or an svg results in an embedded image of the same apparent size.
If you require an exact size in `pt`, consider setting `pt_per_unit = 1`.

Here's an example:

```julia
fig = Figure(resolution = (800, 600))

save("normal.pdf", fig) # size = 600 x 450 pt
save("larger.pdf", fig, pt_per_unit = 2) # size = 1600 x 1200 pt
save("smaller.pdf", fig, pt_per_unit = 0.5) # size = 400 x 300 pt

save("normal.png", fig) # size = 800 x 600 px
save("larger.png", fig, px_per_unit = 2) # size = 1600 x 1200 px
save("smaller.png", fig, px_per_unit = 0.5) # size = 400 x 300 px
```

#### Z-Order

CairoMakie as a 2D engine has no concept of z-clipping, therefore its 3D capabilities are quite limited.
The z-values of 3D plots will have no effect and will be projected flat onto the canvas.
Z-layering is approximated by sorting all plot objects by their z translation value before drawing, after that by parent scene and then insertion order.
Therefore, if you want to draw something on top of something else, but it ends up below, try translating it forward via `translate!(obj, 0, 0, some_positive_z_value)`.
# WGLMakie

[WGLMakie](https://github.com/JuliaPlots/Makie.jl/tree/master/WGLMakie) is the Web-based backend, and is still experimental (though relatively feature-complete). WGLMakie uses [JSServe](https://github.com/SimonDanisch/JSServe.jl) to generate the HTML and JS for the Makie plots.

You can use JSServe and WGLMakie in Pluto, IJulia, Webpages and Documenter to create interactive apps and dashboards, serve them on live webpages, or export them to static HTML.

This tutorial will run through the different modes and what kind of limitations to expect.

First, one should use the new Page mode for anything that displays multiple outputs, like Pluto/IJulia/Documenter.
This creates a single entry point, to connect to the Julia process and load dependencies.
For Documenter, the page needs to be set to `exportable=true, offline=true`.
Exportable has the effect of inlining all data & js dependencies, so that everything can be loaded in a single HTML object.
`offline=true` will make the Page not even try to connect to a running Julia
process, which makes sense for the kind of static export we do in Documenter.


\begin{showhtml}{}
```julia
using JSServe, Markdown
Page(exportable=true, offline=true)
```
\end{showhtml}

After the page got displayed by the frontend, we can start with creating plots and JSServe Apps:


\begin{showhtml}{}
```julia
using WGLMakie
WGLMakie.activate!()
# Set the default resolution to something that fits the Documenter theme
set_theme!(resolution=(800, 400))
scatter(1:4, color=1:4)
```
\end{showhtml}


As you can see, the output is completely static, because we don't have a running Julia server, as it would be the case with e.g. Pluto.
To make the plot interactive, we will need to write more parts of WGLMakie in JS, which is an ongoing effort.
As you can see, the interactivity already keeps working for 3D:

\begin{showhtml}{}
```julia
N = 60
function xy_data(x, y)
    r = sqrt(x^2 + y^2)
    r == 0.0 ? 1f0 : (sin(r)/r)
end
l = range(-10, stop = 10, length = N)
z = Float32[xy_data(x, y) for x in l, y in l]
surface(
    -1..1, -1..1, z,
    colormap = :Spectral
)
```
\end{showhtml}

There are a couple of ways to keep interacting with Plots in a static export.

## Record a statemap

JSServe allows to record a statemap for all widgets, that satisfy the following interface:

```julia
# must be true to be found inside the DOM
is_widget(x) = true
# Updating the widget isn't dependant on any other state (only thing supported right now)
is_independant(x) = true
# The values a widget can iterate
function value_range end
# updating the widget with a certain value (usually an observable)
function update_value!(x, value) end
```

Currently, only sliders overload the interface:

\begin{showhtml}{}
```julia
using Observables

App() do session::Session
    n = 10
    index_slider = Slider(1:n)
    volume = rand(n, n, n)
    slice = map(index_slider) do idx
        return volume[:, :, idx]
    end
    fig = Figure()
    ax, cplot = contour(fig[1, 1], volume)
    rectplot = linesegments!(ax, Rect(-1, -1, 12, 12), linewidth=2, color=:red)
    on(index_slider) do idx
        translate!(rectplot, 0,0,idx)
    end
    heatmap(fig[1, 2], slice)
    slider = DOM.div("z-index: ", index_slider, index_slider.value)
    return JSServe.record_states(session, DOM.div(slider, fig))
end
```
\end{showhtml}

## Execute Javascript directly

JSServe makes it easy to build whole HTML and JS applications.
You can for example directly register javascript function that get run on change.

\begin{showhtml}{}
```julia
using JSServe: onjs

app = App() do session::Session
    s1 = Slider(1:100)
    slider_val = DOM.p(s1[]) # initialize with current value
    # call the `on_update` function whenever s1.value changes in JS:
    onjs(session, s1.value, js"""function on_update(new_value) {
        //interpolating of DOM nodes and other Julia values work mostly as expected:
        const p_element = $(slider_val)
        p_element.innerText = new_value
    }
    """)

    return DOM.div("slider 1: ", s1, slider_val)
end
```
\end{showhtml}

One can also interpolate plots into JS and update those via JS.
The problem is, that there isn't an amazing interface yet.
The returned object is directly a THREE object, with all plot attributes converted into Javascript types.
The good news is, all attributes should be in either `three_scene.material.uniforms`, or `three_scene.geometry.attributes`.
Going forward, we should create an API in WGLMakie, that makes it as easy as in Julia: `plot.attribute = value`.
But while this isn't in place, logging the the returned object makes it pretty easy to figure out what to do - btw, the JS console + logging is amazing and makes it very easy to play around with the object once logged.

\begin{showhtml}{}
```julia
using JSServe: onjs, evaljs, on_document_load

app = App() do session::Session
    s1 = Slider(1:100)
    slider_val = DOM.p(s1[]) # initialize with current value

    fig, ax, splot = scatter(1:4)

    # With on_document_load one can run JS after everything got loaded.
    # This is an alternative to `evaljs`, which we can't use here,
    # since it gets run asap, which means the plots won't be found yet.

    on_document_load(session, js"""
        const plots = $(splot)
        const scatter_plot = plots[0]
        // open the console with ctr+shift+i, to inspect the values
        // tip - you can right click on the log and store the actual variable as a global, and directly interact with it to change the plot.
        console.log(scatter_plot)
        console.log(scatter_plot.material.uniforms)
        console.log(scatter_plot.geometry.attributes)
    """)

    # with the above, we can find out that the positions are stored in `offset`
    # (*sigh*, this is because threejs special cases `position` attributes so it can't be used)
    # Now, lets go and change them when using the slider :)
    onjs(session, s1.value, js"""function on_update(new_value) {
        const plots = $(splot)
        const scatter_plot = plots[0]

        // change first point x + y value
        scatter_plot.geometry.attributes.offset.array[0] = (new_value/100) * 4
        scatter_plot.geometry.attributes.offset.array[1] = (new_value/100) * 4
        // this always needs to be set of geometry attributes after an update
        scatter_plot.geometry.attributes.offset.needsUpdate = true
    }
    """)
    # and for got measures, add a slider to change the color:
    color_slider = Slider(LinRange(0, 1, 100))
    onjs(session, color_slider.value, js"""function on_update(hue) {
        const plot = $(splot)[0]
        const color = new THREE.Color()
        color.setHSL(hue, 1.0, 0.5)
        plot.material.uniforms.color.value.x = color.r
        plot.material.uniforms.color.value.y = color.g
        plot.material.uniforms.color.value.z = color.b
    }""")

    markersize = Slider(1:100)
    onjs(session, markersize.value, js"""function on_update(size) {
        const plot = $(splot)[0]
        plot.material.uniforms.markersize.value.x = size
        plot.material.uniforms.markersize.value.y = size
    }""")
    return DOM.div(s1, color_slider, markersize, fig)
end
```
\end{showhtml}

This summarizes the current state of interactivity with WGLMakie inside static pages.

# Pluto/IJulia

Note that the normal interactivity from Makie is preserved with WGLMakie in e.g. Pluto, as long as the Julia session is running.
Which brings us to setting up Pluto/IJulia sessions! The return value of your first cell must be the return value of the function `Page`.
For example, your first cell can be

```julia
begin
	using JSServe
	Page()
end
```

As is common with files meant to be shared, you might wish to set up a temporary directory so as to not pollute other people's environment. The following code will also be a valid first cell.

```julia
begin
	using Pkg
	Pkg.activate(mktempdir())

	Pkg.add("JSServe")
	using JSServe
	Page()
end
```

If you're accessing the notebook from another PC, you must set:

```julia
begin
	using JSServe
	Page(listen_url="0.0.0.0")
end
```

For more advanced setups consult the `?Page` docs and `JSServe.configure_server!`.

## Styling

You may have noticed, styling isn't really amazing right now.
The good news is, that one can use the whole mighty power of the CSS/HTML universe.
If it wasn't clear so far, JSServe allows to load arbitrary css, and `DOM.xxx` wraps all existing HTML tags.

\begin{showhtml}{}
```julia
using Colors
using JSServe: rows

App() do session::Session

    hue_slider = Slider(0:360)
    color_swatch = DOM.div(class="h-6 w-6 p-2 m-2 rounded shadow")

    onjs(session, hue_slider.value, js"""function (hue){
        $(color_swatch).style.backgroundColor = "hsl(" + hue + ",60%,50%)"
    }""")

    return DOM.div(JSServe.TailwindCSS, rows(hue_slider, color_swatch))
end
```
\end{showhtml}

Tailwind is quite a amazing and has a great documentation especially for CSS beginners:
https://tailwindcss.com/docs/

Note, that JSServe.TailwindCSS is nothing but:

```julia
TailwindCSS = JSServe.Asset("/path/to/tailwind.min.css")
```

So any other CSS file can be used.

It's also pretty easy to make reusable blocks from styled elements.
E.g. the `rows` function above is nothing but:

```julia
rows(args...; class="") = DOM.div(args..., class=class * " flex flex-row")
```

It would be more correct to define it as:

```julia
rows(args...; class="") = DOM.div(JSServe.TailwindCSS, args..., class=class * " flex flex-row")
```

JSServe will then make sure, that `JSServe.TailwindCSS` is loaded, and will only load it once!

Finally, lets create a styled, reusable card componenent:

\begin{showhtml}{}
```julia
using Markdown

struct GridCard
    elements::Any
end

GridCard(elements...) = GridCard(elements)

function JSServe.jsrender(card::GridCard)
    return DOM.div(JSServe.TailwindCSS, card.elements..., class="rounded-lg p-2 m-2 shadow-lg grid auto-cols-max grid-cols-2 gap-4")
end

App() do session::Session
    # We can now use this wherever we want:
    fig = Figure(resolution=(200, 200))
    contour(fig[1,1], rand(4,4))
    card = GridCard(
        Slider(1:100),
        DOM.h1("hello"),
        DOM.img(src="https://julialang.org/assets/infra/logo.svg"),
        fig
    )
    # Markdown creates a DOM as well, and you can interpolate
    # arbitrary jsrender'able elements in there:
    return md"""

    # Wow, Markdown works as well?

    $(card)

    """
end
```
\end{showhtml}

Hopefully, over time there will be helper libraries with lots of stylised elements like the above, to make flashy dashboards with JSServe + WGLMakie.


# Export

Documenter just renders the plots + Page as html,
so if you want to inline WGLMakie/JSServe objects into your own page,
one can just use something like this:

```julia
using WGLMakie, JSServe

using WGLMakie, JSServe
WGLMakie.activate!()

open("index.html", "w") do io
    println(io, """
    <html>
        <head>
        </head>
        <body>
    """)
    # before doing anything else,
    # make sure the Page setup code gets rendered as HTML
    show(io, MIME"text/html"(), Page(exportable=true, offline=true))
    # Then, you can just inline plots or whatever you want :)
    show(io, MIME"text/html"(), scatter(1:4))
    show(io, MIME"text/html"(), surface(rand(4, 4)))
    # or anything else from JSServe, or that can be displayed as html:
    show(io, MIME"text/html"(), JSServe.Slider(1:3))
    println(io, """
        </body>
    </html>
    """)
end
```

# Troubleshooting

## Plots don't display in Safari

Safari users may need to [enable](https://discussions.apple.com/thread/8655829) WebGL.
# GLMakie

[GLMakie](https://github.com/JuliaPlots/Makie.jl/tree/master/GLMakie) is the native, desktop-based backend, and is the most feature-complete.
It requires an OpenGL enabled graphics card with OpenGL version 3.3 or higher.

### Special GLMakie Properties

#### Window Parameters

You can set parameters of the window with the function `set_window_config!` which only takes effect when opening a new window.

```julia
set_window_config!(;
    renderloop = renderloop,
    vsync = false,
    framerate = 30.0,
    float = false,
    pause_rendering = false,
    focus_on_show = false,
    decorated = true,
    title = "Makie"
)
```

## Forcing Dedicated GPU Use In Linux

Normally the dedicated GPU is used for rendering.
If instead an integrated GPU is used, one can tell Julia to use the dedicated GPU while launching julia as `$ sudo DRI_PRIME=1 julia` in the bash terminal.


## Troubleshooting OpenGL

If you get any error loading GLMakie, it likely means, you don't have an OpenGL capable Graphic Card, or you don't have an OpenGL 3.3 capable driver installed.
Note, that most GPUs, even 8 year old integrated ones, support OpenGL 3.3.

On Linux, you can find out your OpenGL version with:
`glxinfo | grep "OpenGL version"`

If you're using an AMD or Intel gpu on linux, you may run into [GLFW#198](https://github.com/JuliaGL/GLFW.jl/issues/198).

If you're on a headless server, you still need to install x-server and
proper graphics drivers.

You can find a demo on how to set that up in this [nextjournal article](https://nextjournal.com/sdanisch/GLMakie-nogpu).

GLMakie's CI has no GPU, so you can also look at [.github/workflows/glmakie.yaml](https://github.com/JuliaPlots/Makie.jl/blob/master/.github/workflows/glmakie.yaml) for a working setup.

If none of these work for you, take a look at the other [backends](/documentation/backends/), which all work without a GPU.

If you get an error pointing to [GLFW.jl](https://github.com/JuliaGL/GLFW.jl), please look into the existing [GLFW issues](https://github.com/JuliaGL/GLFW.jl/issues), and also google for those errors. This is then very likely something that needs fixing in the  [glfw c library](https://github.com/glfw/glfw) or in the GPU drivers.


## WSL setup or X-forwarding

From: [Microsoft/WSL/issues/2855](https://github.com/Microsoft/WSL/issues/2855#issuecomment-358861903)

WSL runs OpenGL alright, but it is not a supported scenario.
From a clean Ubuntu install from the store do:

```
sudo apt install ubuntu-desktop mesa-utils
export DISPLAY=localhost:0
glxgears
```

On the Windows side:

1) install [VcXsrv](https://sourceforge.net/projects/vcxsrv/)
2) choose multiple windows -> display 0 -> start no client -> disable native opengl

Troubleshooting:

1.)  install: `sudo apt-get install -y xorg-dev mesa-utils xvfb libgl1 freeglut3-dev libxrandr-dev libxinerama-dev libxcursor-dev libxi-dev libxext-dev`

2.) WSL has some problems with passing through localhost, so one may need to use: `export DISPLAY=192.168.178.31:0`, with the local ip of the pcs network adapter, which runs VcXsrv

3.) One may need `mv /opt/julia-1.5.2/lib/julia/libstdc++.so.6 /opt/julia-1.5.2/lib/julia/libcpp.backup`, another form of [GLFW#198](https://github.com/JuliaGL/GLFW.jl/issues/198)
# RPRMakie

Experimental ray tracing backend using AMDs [RadeonProRender](https://radeon-pro.github.io/RadeonProRenderDocs/en/index.html).
While it's created by AMD and tailored to Radeon GPUs, it still works just as well for NVidia and Intel GPUs using OpenCL.
It also works on the CPU and even has a hybrid modus to use GPUs and CPUs in tandem to render images.

RadeonProRender currently only works on Windows and Linux, and crashes on OSX when creating the most basic context. If you're on OSX and good at debugging segfaults, help us by debugging:

```julia
using RadeonProRender
RadeonProRender.Context()
```

## Activating and working with RPMakie

To use the backend, just call `activate!` like with all other backends. There are a few extra parameters for RPRMakie:

```julia:docs
# hideall
using RPRMakie, Markdown
println("~~~")
println(Markdown.html(@doc RPRMakie.activate!))
println("~~~")
```
\textoutput{docs}

Since RPRMakie is quite the unique backend and still experimental, there are several gotchas when working with it.

```julia
fig = Figure()
radiance = 10000
# Lights are much more important for ray tracing,
# so most examples will use extra lights and environment lights.
# Note, that RPRMakie is the only backend
# supporting multiple light sources and EnvironmentLights right now
lights = [
    EnvironmentLight(0.5, load(RPR.assetpath("studio026.exr"))),
    PointLight(Vec3f(0, 0, 20), RGBf(radiance, radiance, radiance))]

# Only LScene is supported right now,
# since the other projections don't map to the pysical acurate Camera in RPR.
ax = LScene(fig[1, 1]; scenekw=(lights=lights,))

# to create materials, one needs access to the RPR context.
# Note, if you create an RPRScreen manually, don't display the scene or fig anymore, since that would create a new RPR context, in which resources from the manually created Context would be invalid. Since RPRs error handling is pretty bad, this usually results in Segfaults.
# See below how to render a picture with a manually created context
screen = RPRScreen(ax.scene; iterations=10, plugin=RPR.Northstar)
matsys = screen.matsys
context = screen.context
# You can use lots of materials from RPR.
# Note, that this API may change in the future to a backend  independent representation
# Or at least one, that doesn't need to access the RPR context
mat = RPR.Chrome(matsys)
# The material attribute is specific to RPRMakie and gets ignored by other Backends. This may change in the future
mesh!(ax, Sphere(Point3f, 0), material=mat)

# There are three main ways to turn a Makie scene into a picture:
# Get the colorbuffer of the RPRScreen. RPRScreen also has `show` overloaded for the mime `image\png` so it should display in IJulia/Jupyter/VSCode.
image = colorbuffer(screen)::Matrix{RGB{N0f8}}
# Replace a specific (sub) LScene with RPR, and display the whole scene interactively in GLMakie
using GLMakie
refres = Observable(nothing) # Optional observable that triggers
GLMakie.activate!(); display(fig) # Make sure to display scene first in GLMakie
# Replace the scene with an interactively rendered RPR output.
# See more about this in the GLMakie interop example
context, task = RPRMakie.replace_scene_rpr!(ax.scene, screen; refresh=refresh)
# If one doesn't create the RPRScreen manually to create custom materials,
# display(ax.scene), show(io, MIME"image/png", ax.scene), save("rpr.png", ax.scene)
# Should work just like with other backends.
# Note, that only the scene from LScene can be displayed directly, but soon, `display(fig)` should also work.
```


There are several examples showing different aspects of how to use RPRMakie.
The examples are in [RPRMakie/examples](https://github.com/JuliaPlots/Makie.jl/tree/master/RPRMakie/examples)

## MaterialX and predefined materials (materials.jl)

There are several predefined materials one can use in RadeonProRender.
RPR also supports the [MaterialX](https://www.materialx.org/) standard to load a wide range of predefined Materials. Make sure to use the Northstar backend for `Matx` materials.

~~~
<input id="hidecode5" class="hidecode" type="checkbox">
~~~
```julia
using GeometryBasics, RPRMakie
using Colors, FileIO
using Colors: N0f8

radiance = 500
lights = [EnvironmentLight(1.0, load(RPR.assetpath("studio026.exr"))),
            PointLight(Vec3f(10), RGBf(radiance, radiance, radiance * 1.1))]
fig = Figure(; resolution=(1500, 700))
ax = LScene(fig[1, 1]; show_axis=false, scenekw=(lights=lights,))
screen = RPRScreen(ax.scene; plugin=RPR.Northstar, iterations=400)

matsys = screen.matsys
emissive = RPR.EmissiveMaterial(matsys)
diffuse = RPR.DiffuseMaterial(matsys)
glass = RPR.Glass(matsys)
plastic = RPR.Plastic(matsys)
chrome = RPR.Chrome(matsys)
dielectric = RPR.DielectricBrdfX(matsys)
gold = RPR.SurfaceGoldX(matsys)

materials = [glass chrome;
                gold dielectric;
                emissive plastic]

mesh!(ax, load(Makie.assetpath("matball_floor.obj")); color=:white)
palette = reshape(Makie.default_palettes.color[][1:6], size(materials))

for i in CartesianIndices(materials)
    x, y = Tuple(i)
    mat = materials[i]
    mplot = if mat === emissive
        matball!(ax, diffuse; inner=emissive, color=nothing)
    else
        matball!(ax, mat; color=nothing)
    end
    v = Vec3f(((x, y) .- (0.5 .* size(materials)) .- 0.5)..., 0)
    translate!(mplot, 0.9 .* (v .- Vec3f(0, 3, 0)))
end
cam = cameracontrols(ax.scene)
cam.eyeposition[] = Vec3f(-0.3, -5.5, 0.9)
cam.lookat[] = Vec3f(0.5, 0, -0.5)
cam.upvector[] = Vec3f(0, 0, 1)
cam.fov[] = 35
emissive.color = Vec3f(4, 2, 2)
image = colorbuffer(screen)
save("materials.png", image)
```
~~~
<label for="hidecode5" class="hidecode"></label>
~~~

~~~
<img src="/assets/materials.png">
~~~

## Advanced custom material (earth_topography.jl)

~~~
<input id="hidecode4" class="hidecode" type="checkbox">
~~~
```julia
using NCDatasets, ColorSchemes, RPRMakie
using ImageShow, FileIO

# Taken from https://lazarusa.github.io/BeautifulMakie/GeoPlots/topography/
cmap = dataset = Dataset(joinpath(@__DIR__, "ETOPO1_halfdegree.nc"))
lon = dataset["lon"][:]
lat = dataset["lat"][:]
data = Float32.(dataset["ETOPO1avg"][:, :])

function glow_material(data_normed)
    emission_weight = map(data_normed) do i
        return Float32(i < 0.7 ? 0.0 : i)
    end
    emission_color = map(data_normed) do i
        em = i * 2
        return RGBf(em * 2.0, em * 0.4, em * 0.3)
    end

    return (
        reflection_weight = 1,
        reflection_color = RGBf(0.5, 0.5, 1.0),
        reflection_metalness = 0,
        reflection_ior = 1.4,
        diffuse_weight = 1,
        emission_weight = emission_weight',
        emission_color = emission_color',
    )
end

RPRMakie.activate!(iterations=32, plugin=RPR.Northstar)
fig = Figure(; resolution=(2000, 800))
radiance = 30000
lights = [EnvironmentLight(1.0, load(RPR.assetpath("studio026.exr"))),
            PointLight(Vec3f(0, 100, 100), RGBf(radiance, radiance, radiance))]

ax = LScene(fig[1, 1]; show_axis=false, scenekw=(lights=lights,))

mini, maxi = extrema(data)
data_normed = ((data .- mini) ./ (maxi - mini))

material = glow_material(data_normed)

pltobj = surface!(ax, lon, lat, data_normed .* 20;
                    material=material, colormap=[:black, :white, :brown],
                    colorrange=(0.2, 0.8) .* 20)
# Set the camera to a nice angle
cam = cameracontrols(ax.scene)
cam.eyeposition[] = Vec3f(3, -300, 300)
cam.lookat[] = Vec3f(0)
cam.upvector[] = Vec3f(0, 0, 1)
cam.fov[] = 23

save("topographie.png", ax.scene)
```
~~~
<label for="hidecode4" class="hidecode"></label>
~~~

~~~
<img src="/assets/topographie.png">
~~~

## GLMakie interop (opengl_interop.jl)

RPRMakie doesn't support layouting and sub scenes yet, but you can replace a single scene with a RPR renedered, interactive window.
This is especially handy, to show 2d graphics and interactive UI elements next to a ray traced scene and interactively tune camera and material parameters.

~~~
<input id="hidecode3" class="hidecode" type="checkbox">
~~~
```julia
using GLMakie, GeometryBasics, RPRMakie, RadeonProRender
using Colors, FileIO
using Colors: N0f8

f = (u, v) -> cos(v) * (6 - (5 / 4 + sin(3 * u)) * sin(u - 3 * v))
g = (u, v) -> sin(v) * (6 - (5 / 4 + sin(3 * u)) * sin(u - 3 * v))
h = (u, v) -> -cos(u - 3 * v) * (5 / 4 + sin(3 * u));
u = range(0; stop=2π, length=150)
v = range(0; stop=2π, length=150)
radiance = 500
lights = [EnvironmentLight(1.0, load(RPR.assetpath("studio026.exr"))),
          PointLight(Vec3f(10), RGBf(radiance, radiance, radiance * 1.1))]

fig = Figure(; resolution=(1500, 1000))
ax = LScene(fig[1, 1]; show_axis=false, scenekw=(lights=lights,))
screen = RPRMakie.RPRScreen(size(ax.scene); plugin=RPR.Tahoe)
material = RPR.UberMaterial(screen.matsys)

surface!(ax, f.(u, v'), g.(u, v'), h.(u, v'); ambient=Vec3f(0.5), diffuse=Vec3f(1), specular=0.5,
         colormap=:balance, material=material)

function Input(fig, val::RGB)
    hue = Slider(fig; range=1:380, width=200)
    lightness = Slider(fig; range=LinRange(0, 1, 100), width=200)
    labels = [Label(fig, "hue"; halign=:left), Label(fig, "light"; halign=:left)]
    layout = grid!(hcat(labels, [hue, lightness]))
    hsl = HSL(val)
    set_close_to!(hue, hsl.h)
    set_close_to!(lightness, hsl.l)
    color = map((h, l) -> RGB(HSL(h, 0.9, l)), hue.value, lightness.value)
    return color, layout
end

function Input(fig, val::Vec4)
    s = Slider(fig; range=LinRange(0, 1, 100), width=200)
    set_close_to!(s, first(val))
    return map(x -> Vec4f(x), s.value), s
end

function Input(fig, val::Bool)
    toggle = Toggle(fig; active=val)
    return toggle.active, toggle
end

sliders = (reflection_color=Input(fig, RGB(0, 0, 0)), reflection_weight=Input(fig, Vec4(0)),
           reflection_roughness=Input(fig, Vec4(0)), reflection_anisotropy=Input(fig, Vec4(0)),
           reflection_anisotropy_rotation=Input(fig, Vec4(0)), reflection_mode=Input(fig, Vec4(0)),
           reflection_ior=Input(fig, Vec4(0)), reflection_metalness=Input(fig, Vec4(0)),
           refraction_color=Input(fig, RGB(0, 0, 0)), refraction_weight=Input(fig, Vec4(0)),
           refraction_roughness=Input(fig, Vec4(0)), refraction_ior=Input(fig, Vec4(0)),
           refraction_absorption_color=Input(fig, RGB(0, 0, 0)),
           refraction_absorption_distance=Input(fig, Vec4(0)), refraction_caustics=Input(fig, true),
           sss_scatter_color=Input(fig, RGB(0, 0, 0)), sss_scatter_distance=Input(fig, Vec4(0)),
           sss_scatter_direction=Input(fig, Vec4(0)), sss_weight=Input(fig, Vec4(0)),
           sss_multiscatter=Input(fig, false), backscatter_weight=Input(fig, Vec4(0)),
           backscatter_color=Input(fig, RGB(0, 0, 0)))

labels = []
inputs = []
refresh = Observable(nothing)
for (key, (obs, input)) in pairs(sliders)
    push!(labels, Label(fig, string(key); align=:left))
    push!(inputs, input)
    on(obs) do value
        @show key value
        setproperty!(material, key, value)
        return notify(refresh)
    end
end

fig[1, 2] = grid!(hcat(labels, inputs); width=500)
GLMakie.activate!()

cam = cameracontrols(ax.scene)
cam.eyeposition[] = Vec3f(22, 0, 17)
cam.lookat[] = Vec3f(0, 0, -1)
cam.upvector[] = Vec3f(0, 0, 1)
cam.fov[] = 30

display(fig)

context, task = RPRMakie.replace_scene_rpr!(ax.scene, screen; refresh=refresh)

# Change light parameters interactively
begin
    lights[1].intensity[] = 1.5
    lights[2].radiance[] = RGBf(1000, 1000, 1000)
    lights[2].position[] = Vec3f(3, 10, 10)
    notify(refresh)
end
```
~~~
<label for="hidecode3" class="hidecode"></label>
~~~

~~~
<video autoplay controls src="/assets/opengl_interop.mp4">
</video>
~~~

## Animations (lego.jl)

Not all objects support updating via Observables yet, but translations, camera etc are already covered and can be used together with Makie's standard animation API.

~~~
<input id="hidecode2" class="hidecode" type="checkbox">
~~~
```julia
# Example inspiration and Lego model by https://github.com/Kevin-Mattheus-Moerman
# https://twitter.com/KMMoerman/status/1417759722963415041
using MeshIO, FileIO, GeometryBasics, RPRMakie

colors = Dict(
    "eyes" => "#000",
    "belt" => "#000059",
    "arm" => "#009925",
    "leg" => "#3369E8",
    "torso" => "#D50F25",
    "head" => "yellow",
    "hand" => "yellow"
)

origins = Dict(
    "arm_right" => Point3f(0.1427, -6.2127, 5.7342),
    "arm_left" => Point3f(0.1427, 6.2127, 5.7342),
    "leg_right" => Point3f(0, -1, -8.2),
    "leg_left" => Point3f(0, 1, -8.2),
)

rotation_axes = Dict(
    "arm_right" => Vec3f(0.0000, -0.9828, 0.1848),
    "arm_left" => Vec3f(0.0000, 0.9828, 0.1848),
    "leg_right" => Vec3f(0, -1, 0),
    "leg_left" => Vec3f(0, 1, 0),
)

function plot_part!(scene, parent, name::String)
    m = load(assetpath("lego_figure_" * name * ".stl"))
    color = colors[split(name, "_")[1]]
    trans = Transformation(parent)
    ptrans = Makie.transformation(parent)
    origin = get(origins, name, nothing)
    if !isnothing(origin)
        centered = m.position .- origin
        m = GeometryBasics.Mesh(meta(centered; normals=m.normals), faces(m))
        translate!(trans, origin)
    else
        translate!(trans, -ptrans.translation[])
    end
    return mesh!(scene, m; color=color, transformation=trans)
end

function plot_lego_figure(s, floor=true)
    # Plot hierarchical mesh!
    figure = Dict()
    # Plot hierarchical mesh!
    figure["torso"] = plot_part!(s, s, "torso")
        figure["head"] = plot_part!(s, figure["torso"], "head")
            figure["eyes_mouth"] = plot_part!(s, figure["head"], "eyes_mouth")
        figure["arm_right"] = plot_part!(s, figure["torso"], "arm_right")
            figure["hand_right"] = plot_part!(s, figure["arm_right"], "hand_right")
        figure["arm_left"] = plot_part!(s, figure["torso"], "arm_left")
            figure["hand_left"] = plot_part!(s, figure["arm_left"], "hand_left")
        figure["belt"] = plot_part!(s, figure["torso"], "belt")
            figure["leg_right"] = plot_part!(s, figure["belt"], "leg_right")
            figure["leg_left"] = plot_part!(s, figure["belt"], "leg_left")
    # lift the little guy up
    translate!(figure["torso"], 0, 0, 20)
    # add some floor
    floor && mesh!(s, Rect3f(Vec3f(-400, -400, -2), Vec3f(800, 800, 2)), color=:white)
    return figure
end

RPRMakie.activate!(iterations=200, plugin=RPR.Northstar)
radiance = 50000
lights = [
    EnvironmentLight(1.5, rotl90(load(assetpath("sunflowers_1k.hdr"))')),
    PointLight(Vec3f(50, 0, 200), RGBf(radiance, radiance, radiance*1.1)),
]
s = Scene(resolution=(500, 500), lights=lights)

cam3d!(s)
c = cameracontrols(s)
c.near[] = 5
c.far[] = 1000
update_cam!(s, c, Vec3f(100, 30, 80), Vec3f(0, 0, -10))
figure = plot_lego_figure(s)

rot_joints_by = 0.25*pi
total_translation = 50
animation_strides = 10

a1 = LinRange(0, rot_joints_by, animation_strides)
angles = [a1; reverse(a1[1:end-1]); -a1[2:end]; reverse(-a1[1:end-1]);]
nsteps = length(angles); #Number of animation steps
translations = LinRange(0, total_translation, nsteps)

Makie.record(s, "lego_walk.mp4", zip(translations, angles)) do (translation, angle)

    # Rotate right arm + hand
    for name in ["arm_left", "arm_right", "leg_left", "leg_right"]
        rotate!(figure[name], rotation_axes[name], angle)
    end
    translate!(figure["torso"], translation, 0, 20)
end
```
~~~
<label for="hidecode2" class="hidecode"></label>
~~~

~~~
<video autoplay controls src="/assets/lego_walk.mp4">
</video>
~~~

## Earth example

~~~
<input id="hidecode1" class="hidecode" type="checkbox">
~~~
```julia
# by Lazaro Alonso
# taken from: https://lazarusa.github.io/BeautifulMakie/GeoPlots/submarineCables3D/
using GeoMakie, Downloads
using GeoJSON, GeoInterface
using FileIO
using RPRMakie
# data from
# https://github.com/telegeography/www.submarinecablemap.com
urlPoints = "https://raw.githubusercontent.com/telegeography/www.submarinecablemap.com/master/web/public/api/v3/landing-point/landing-point-geo.json"
urlCables = "https://raw.githubusercontent.com/telegeography/www.submarinecablemap.com/master/web/public/api/v3/cable/cable-geo.json"

landPoints = Downloads.download(urlPoints, IOBuffer())
landCables = Downloads.download(urlCables, IOBuffer())

land_geoPoints = GeoJSON.read(seekstart(landPoints))
land_geoCables = GeoJSON.read(seekstart(landCables))

toPoints = GeoMakie.geo2basic(land_geoPoints)
feat = GeoInterface.features(land_geoCables)
toLines = GeoInterface.coordinates.(GeoInterface.geometry.(feat))

# broken lines at -180 and 180... they should
# be the same line and be in the same array.

# some 3D transformations
function toCartesian(lon, lat; r = 1.02, cxyz = (0, 0, 0))
    x = cxyz[1] + r * cosd(lat) * cosd(lon)
    y = cxyz[2] + r * cosd(lat) * sind(lon)
    z = cxyz[3] + r * sind(lat)
    return (x, y, z)
end

toPoints3D = [Point3f([toCartesian(point[1], point[2])...]) for point in toPoints]

splitLines3D = []
for i in 1:length(toLines)
    for j in 1:length(toLines[i])
        ptsLines = toLines[i][j]
        tmp3D = []
        for k in 1:length(ptsLines)
            x, y = ptsLines[k]
            x, y, z = toCartesian(x, y)
            push!(tmp3D, [x, y, z])
        end
        push!(splitLines3D, Point3f.(tmp3D))
    end
end

earth_img = load(Downloads.download("https://upload.wikimedia.org/wikipedia/commons/5/56/Blue_Marble_Next_Generation_%2B_topography_%2B_bathymetry.jpg"))
Makie.inline!(true)
# the actual plot !
RPRMakie.activate!(; iterations=100)
scene = with_theme(theme_dark()) do
    fig = Figure(; resolution=(1000, 1000))
    radiance = 30
    lights = [EnvironmentLight(0.5, load(RPR.assetpath("starmap_4k.tif"))),
              PointLight(Vec3f(1, 1, 3), RGBf(radiance, radiance, radiance))]
    ax = LScene(fig[1, 1]; show_axis=false, scenekw=(lights=lights,))
    n = 1024 ÷ 4 # 2048
    θ = LinRange(0, pi, n)
    φ = LinRange(-pi, pi, 2 * n)
    xe = [cos(φ) * sin(θ) for θ in θ, φ in φ]
    ye = [sin(φ) * sin(θ) for θ in θ, φ in φ]
    ze = [cos(θ) for θ in θ, φ in φ]
    surface!(ax, xe, ye, ze; color=earth_img)
    meshscatter!(toPoints3D; color=1:length(toPoints3D), markersize=0.005, colormap=:plasma)
    colors = Makie.default_palettes.color[]
    c = Iterators.cycle(colors)
    foreach(((l, c),) -> lines!(ax, l; linewidth=2, color=c), zip(splitLines3D, c))
    ax.scene.camera_controls.eyeposition[] = Vec3f(1.5)
    return ax.scene
end

save("submarin_cables.png", scene)
```
~~~
<label for="hidecode1" class="hidecode"></label>
~~~


~~~
<img src="/assets/earth.png">
~~~
# Layoutables & Widgets

All examples in this section are presented as static CairoMakie vector graphics for clarity of visuals.
Keep in mind that CairoMakie is not interactive.
Use GLMakie for interactive widgets, as WGLMakie currently doesn't have picking implemented.

## Example gallery

{{list_folder_with_images layoutables}}# Plotting functions

{{list_folder_with_images plotting_functions}}# IntervalSlider

The interval slider selects an interval (low, high) from the supplied attribute `range`.
The (approximate) start values can be set with `startvalues`.

The currently selected interval is in the attribute `interval` and is a Tuple of `(low, high)`.
Don't change this value manually, but use the function `set_close_to!(intslider, v1, v2)`.
This is necessary to ensure the values are actually present in the `range` attribute.

You can click anywhere outside of the currently selected range and the closer interval edge will jump to the point.
You can then drag the edge around.
When hovering over the slider, the larger button indicates the edge that will react.

If the mouse hovers over the central area of the interval and both buttons are enlarged, clicking and dragging shifts the interval around as a whole.

You can double-click the slider to reset it to the values present in `startvalues`.
If `startvalues === Makie.automatic`, the full interval will be selected (this is the default).

If you set the attribute `snap = false`, the slider will move continously while dragging and only jump to the closest available values when releasing the mouse.

\begin{examplefigure}{}
```julia
using CairoMakie
Makie.inline!(true) # hide
CairoMakie.activate!() # hide

f = Figure()
Axis(f[1, 1], limits = (0, 1, 0, 1))

rs_h = IntervalSlider(f[2, 1], range = LinRange(0, 1, 1000),
    startvalues = (0.2, 0.8))
rs_v = IntervalSlider(f[1, 2], range = LinRange(0, 1, 1000),
    startvalues = (0.4, 0.9), horizontal = false)

labeltext1 = lift(rs_h.interval) do int
    string(round.(int, digits = 2))
end
Label(f[3, 1], labeltext1, tellwidth = false)
labeltext2 = lift(rs_v.interval) do int
    string(round.(int, digits = 2))
end
Label(f[1, 3], labeltext2,
    tellheight = false, rotation = pi/2)

points = rand(Point2f, 300)

# color points differently if they are within the two intervals
colors = lift(rs_h.interval, rs_v.interval) do h_int, v_int
    map(points) do p
        (h_int[1] < p[1] < h_int[2]) && (v_int[1] < p[2] < v_int[2])
    end
end

scatter!(points, color = colors, colormap = [:black, :orange], strokewidth = 0)

f
```
\end{examplefigure}

# Slider

A simple slider without a label. You can create a label using a `Label` object,
for example. You need to specify a range that constrains the slider's possible values.

The currently selected value is in the attribute `value`.
Don't change this value manually, but use the function `set_close_to!(slider, value)`.
This is necessary to ensure the value is actually present in the `range` attribute.

You can double-click the slider to reset it (approximately) to the value present in `startvalue`.

If you set the attribute `snap = false`, the slider will move continously while dragging and only jump to the closest available value when releasing the mouse.

\begin{examplefigure}{}
```julia
using GLMakie
GLMakie.activate!() # hide
fig = Figure()

ax = Axis(fig[1, 1])

sl_x = Slider(fig[2, 1], range = 0:0.01:10, startvalue = 3)
sl_y = Slider(fig[1, 2], range = 0:0.01:10, horizontal = false, startvalue = 6)

point = lift(sl_x.value, sl_y.value) do x, y
    Point2f(x, y)
end

scatter!(point, color = :red, markersize = 20)

limits!(ax, 0, 10, 0, 10)

fig
```
\end{examplefigure}


## Labelled sliders and grids

To create a horizontal layout containing a label, a slider, and a value label, use the convenience function \apilink{labelslider!}.

If you need multiple aligned rows of sliders, use \apilink{labelslidergrid!}.

The column with the value labels is automatically set to a fixed width, so that the layout doesn't jitter when sliders are dragged and the value labels change their widths.
This width is chosen by setting each slider to a few values and recording the maximum label width.
Alternatively, you can set the width manually with the keyword argument `value_column_width`.

You can pass either format functions or format strings as used by `Formatting.format`.

\begin{examplefigure}{}
```julia
using GLMakie

fig = Figure()

ax = Axis(fig[1, 1])

lsgrid = labelslidergrid!(
    fig,
    ["Voltage", "Current", "Resistance"],
    [0:0.1:10, 0:0.1:20, 0:0.1:30];
    formats = "{:.1f}" .* ["V", "A", "Ω"],
    width = 350,
    tellheight = false)

fig[1, 2] = lsgrid.layout

sliderobservables = [s.value for s in lsgrid.sliders]
bars = lift(sliderobservables...) do slvalues...
    [slvalues...]
end

barplot!(ax, bars, color = [:yellow, :orange, :red])
ylims!(ax, 0, 30)

set_close_to!(lsgrid.sliders[1], 5.3)
set_close_to!(lsgrid.sliders[2], 10.2)
set_close_to!(lsgrid.sliders[3], 15.9)

fig
```
\end{examplefigure}


# Menu

A dropdown menu with `options`, where each element's label is determined with `optionlabel(element)`
and the value with `optionvalue(element)`. The default behavior is to treat a 2-element tuple
as `(label, value)` and any other object as `value`, where `label = string(value)`.

The attribute `selection` is set to `optionvalue(element)` when the element's entry is selected.

\begin{examplefigure}{}
```julia
using GLMakie
GLMakie.activate!() # hide
fig = Figure()

menu = Menu(fig, options = ["viridis", "heat", "blues"])

funcs = [sqrt, x->x^2, sin, cos]

menu2 = Menu(fig, options = zip(["Square Root", "Square", "Sine", "Cosine"], funcs))

fig[1, 1] = vgrid!(
    Label(fig, "Colormap", width = nothing),
    menu,
    Label(fig, "Function", width = nothing),
    menu2;
    tellheight = false, width = 200)

ax = Axis(fig[1, 2])

func = Observable{Any}(funcs[1])

ys = lift(func) do f
    f.(0:0.3:10)
end
scat = scatter!(ax, ys, markersize = 10px, color = ys)

cb = Colorbar(fig[1, 3], scat)

on(menu.selection) do s
    scat.colormap = s
end

on(menu2.selection) do s
    func[] = s
    autolimits!(ax)
end

menu2.is_open = true

fig
```
\end{examplefigure}

## Menu direction

You can change the direction of the menu with `direction = :up` or `direction = :down`. By default, the direction is determined automatically to avoid cutoff at the figure boundaries.


\begin{examplefigure}{}
```julia
using GLMakie
GLMakie.activate!() # hide

fig = Figure()

menu = Menu(fig[1, 1], options = ["A", "B", "C"])
menu2 = Menu(fig[3, 1], options = ["A", "B", "C"])

menu.is_open = true
menu2.is_open = true

fig
```
\end{examplefigure}


# Legend

## Creating A Legend From Elements

You can create a basic Legend by passing a vector of legend entries and a vector of labels, plus an optional title as the third argument.

The elements in the vector of legend entries can either be plot objects or LegendElements like LineElement, MarkerElement and PolyElement.
Or they can be vectors of such objects that will be layered together as one.

### Legend element attributes

The standard plot objects like `Scatter` or `Lines` have predefined conversions to `MarkerElement`s and `LineElement`s that copy the relevant plot attributes to the legend element.
If an attribute has a vector-like value, it falls back to the scalar default of the legend.
The legend defaults themselves are by default inherited from the main theme.
For example, `polystrokewidth` of the legend falls back to `patchstrokewidth` of the main theme.
In the following example, you can see that the legend for `sca2` copies the `:rect` marker but not the vector-valued color.

\begin{examplefigure}{svg = true}
```julia
using CairoMakie
CairoMakie.activate!() # hide

f = Figure()

Axis(f[1, 1])

xs = 0:0.5:10
ys = sin.(xs)
lin = lines!(xs, ys, color = :blue)
sca = scatter!(xs, ys, color = :red)
sca2 = scatter!(xs, ys .+ 0.5, color = 1:length(xs), marker = :rect)

Legend(f[1, 2],
    [lin, sca, [lin, sca], sca2],
    ["a line", "some dots", "both together", "rect markers"])

f
```
\end{examplefigure}
## Creating A Legend From An Axis

You can also create a Legend by passing it an axis object, like `Axis`, `LScene` or `Scene`.
All plots that have a `label` attribute set will be put into the legend, in the order that they appear in the axis, and you can optionally pass a title as the third argument.

\begin{examplefigure}{svg = true}
```julia
using CairoMakie

f = Figure()

ax = f[1, 1] = Axis(f)

lines!(0..15, sin, label = "sin", color = :blue)
lines!(0..15, cos, label = "cos", color = :red)
lines!(0..15, x -> -cos(x), label = "-cos", color = :green)

f[1, 2] = Legend(f, ax, "Trig Functions", framevisible = false)

f
```
\end{examplefigure}
With the keywords `merge` and `unique` you can control how plot objects with the same labels are treated.
If `merge` is `true`, all plot objects with the same label will be layered on top of each other into one legend entry.
If `unique` is `true`, all plot objects with the same plot type and label will be reduced to one occurrence.

\begin{examplefigure}{svg = true}
```julia
using CairoMakie

f = Figure()

traces = cumsum(randn(10, 5), dims = 1)

for (i, (merge, unique)) in enumerate(
        Iterators.product([false, true], [false true]))

    axis = Axis(f[fldmod1(i, 2)...],
        title = "merge = $merge, unique = $unique")

    for trace in eachcol(traces)
        lines!(trace, label = "single", color = (:black, 0.2))
    end

    mu = vec(sum(traces, dims = 2) ./ 5)
    lines!(mu, label = "mean")
    scatter!(mu, label = "mean")

    axislegend(axis, merge = merge, unique = unique)

end

f
```
\end{examplefigure}
## Multi-Bank Legend

You can control the number of banks with the `nbanks` attribute. Banks are columns
when in vertical mode, and rows when in horizontal mode.

\begin{examplefigure}{svg = true}
```julia
using CairoMakie

f = Figure()

Axis(f[1, 1])

xs = 0:0.1:10
lins = [lines!(xs, sin.(xs .+ 3v), color = RGBf(v, 0, 1-v)) for v in 0:0.1:1]

Legend(f[1, 2], lins, string.(1:length(lins)), nbanks = 3)

f
```
\end{examplefigure}
## Legend Inside An Axis

The `axislegend` function is a quick way to add a legend to an Axis.
You can pass a selected axis plus arguments which are forwarded to the `Legend` constructor, or the current axis is used by default.
If you pass only a string, it's used as the title with the current axis.

The position can be set via a shortcut symbol, first halign (l, r, c) then valign (b, t, c), such as :lt for left, top and :cb for center bottom.

\begin{examplefigure}{svg = true}
```julia
using CairoMakie

f = Figure()

ax = Axis(f[1, 1])

sc1 = scatter!(randn(10, 2), color = :red, label = "Red Dots")
sc2 = scatter!(randn(10, 2), color = :blue, label = "Blue Dots")
scatter!(randn(10, 2), color = :orange, label = "Orange Dots")
scatter!(randn(10, 2), color = :cyan, label = "Cyan Dots")

axislegend()

axislegend("Titled Legend", position = :lb)

axislegend(ax, [sc1, sc2], ["One", "Two"], "Selected Dots", position = :rb,
    orientation = :horizontal)

f
```
\end{examplefigure}
Alternatively, you can simply add a Legend to the same layout slot
that an axis lives in. As long as the axis is bigger than the legend you can
set the legend's `tellheight` and `tellwidth` to `false` and position it using the align
variables. You can use the margin keyword to keep the legend from touching the axis
spines.

\begin{examplefigure}{svg = true}
```julia
using CairoMakie

haligns = [:left, :right, :center]
valigns = [:top, :bottom, :center]

f = Figure()

Axis(f[1, 1])

xs = 0:0.1:10
lins = [lines!(xs, sin.(xs .* i), color = color)
    for (i, color) in zip(1:3, [:red, :blue, :green])]

for (j, ha, va) in zip(1:3, haligns, valigns)
    Legend(
        f[1, 1], lins, ["Line $i" for i in 1:3],
        "$ha & $va",
        tellheight = false,
        tellwidth = false,
        margin = (10, 10, 10, 10),
        halign = ha, valign = va, orientation = :horizontal
    )
end

f
```
\end{examplefigure}
## Creating Legend Entries Manually

Sometimes you might want to construct legend entries from scratch to have maximum
control. So far you can use `LineElement`s, `MarkerElement`s or `PolyElement`s.
The attributes for these elements are the following (the `[]` parts can be left out when constructing these elements directly, but have to be fully written out for the attributes that the legend holds):

```julia
# LineElement
[line]points, [line]color, linestyle, linewidth

# MarkerElement
[marker]points, marker, markersize, [marker]color,
[marker]strokewidth, [marker]strokecolor

# PolyElement
[poly]points, [poly]color, [poly]strokewidth, [poly]strokecolor
```

The attributes `linepoints`, `markerpoints` and `polypoints` decide where in the legend entry patch rectangle the plot objects are placed.
These values should be normalized to a 1 by 1 rectangle, and the final shape depends on the `patchsize` of the legend.
For example, if you want wider line and poly markers, you could set the `patchsize` of the legend to `(50, 30)`.

\begin{examplefigure}{svg = true}
```julia
using CairoMakie

f = Figure()

Axis(f[1, 1])

elem_1 = [LineElement(color = :red, linestyle = nothing),
          MarkerElement(color = :blue, marker = 'x', markersize = 15,
          strokecolor = :black)]

elem_2 = [PolyElement(color = :red, strokecolor = :blue, strokewidth = 1),
          LineElement(color = :black, linestyle = :dash)]

elem_3 = LineElement(color = :green, linestyle = nothing,
        points = Point2f[(0, 0), (0, 1), (1, 0), (1, 1)])

elem_4 = MarkerElement(color = :blue, marker = 'π', markersize = 15,
        points = Point2f[(0.2, 0.2), (0.5, 0.8), (0.8, 0.2)])

elem_5 = PolyElement(color = :green, strokecolor = :black, strokewidth = 2,
        points = Point2f[(0, 0), (1, 0), (0, 1)])

Legend(f[1, 2],
    [elem_1, elem_2, elem_3, elem_4, elem_5],
    ["Line & Marker", "Poly & Line", "Line", "Marker", "Poly"],
    patchsize = (35, 35), rowgap = 10)

f
```
\end{examplefigure}
## Horizontal Legend

In case you want the legend entries to be listed horizontally, set the `orientation`
attribute to `:horizontal`. In this case the `nbanks` attribute refers to the
number of rows instead of columns. To keep an adjacent axis from potentially shrinking to
the width of the horizontal legend, set `tellwidth = false` and `tellheight = true`
if you place the legend below or above the axis.

\begin{examplefigure}{svg = true}
```julia
using CairoMakie

f = Figure()

Axis(f[1, 1])

xs = 0:0.5:10
ys = sin.(xs)
lin = lines!(xs, ys, color = :blue)
sca = scatter!(xs, ys, color = :red, markersize = 15)

Legend(f[1, 2], [lin, sca, lin], ["a line", "some dots", "line again"])

Legend(f[2, 1], [lin, sca, lin], ["a line", "some dots", "line again"],
    orientation = :horizontal, tellwidth = false, tellheight = true)

f
```
\end{examplefigure}
## Multi-Group Legends

Sometimes a legend consists of multiple groups, for example in a plot where both
marker size and color are varied and those properties need to be visualized
separately, but still together in one legend. Each group's content is given as
an array of elements and an array of labels, each within one collective array.
You can shift the position of the titles relative to each group with the
`titleposition` attribute, either `:left` or `:top`.

\begin{examplefigure}{svg = true}
```julia
using CairoMakie

f = Figure()

Axis(f[1, 1])

markersizes = [5, 10, 15, 20]
colors = [:red, :green, :blue, :orange]

for ms in markersizes, color in colors
    scatter!(randn(5, 2), markersize = ms, color = color)
end

group_size = [MarkerElement(marker = :circle, color = :black,
    strokecolor = :transparent,
    markersize = ms) for ms in markersizes]

group_color = [PolyElement(color = color, strokecolor = :transparent)
    for color in colors]

legends = [Legend(f,
    [group_size, group_color],
    [string.(markersizes), string.(colors)],
    ["Size", "Color"]) for _ in 1:6]

f[1, 2:4] = legends[1:3]
f[2:4, 2] = legends[4:6]

for l in legends[4:6]
    l.orientation = :horizontal
    l.tellheight = true
    l.tellwidth = false
end

legends[2].titleposition = :left
legends[5].titleposition = :left

legends[3].nbanks = 2
legends[5].nbanks = 2
legends[6].nbanks = 2

f
```
\end{examplefigure}

# Textbox

The `Textbox` supports entry of a simple, single-line string, with optional validation logic.

\begin{examplefigure}{}
```julia
using GLMakie
using CairoMakie # hide
CairoMakie.activate!() # hide

f = Figure()
Textbox(f[1, 1], placeholder = "Enter a string...")
Textbox(f[2, 1], width = 300)

f
```
\end{examplefigure}

## Validation

The `validator` attribute is used with `validate_textbox(string, validator)` to determine if the current string is valid. It can be a `Regex` that needs to match the complete string, or a `Function` taking a `String` as input and returning a `Bool`. If the validator is a type T (for example `Float64`), validation will be `tryparse(string, T)`. The textbox will not allow submitting the currently entered value if the validator doesn't pass.

\begin{examplefigure}{}
```julia
using GLMakie
using CairoMakie # hide
CairoMakie.activate!() # hide

f = Figure()

tb = Textbox(f[2, 1], placeholder = "Enter a frequency",
    validator = Float64, tellwidth = false)

frequency = Observable(1.0)

on(tb.stored_string) do s
    frequency[] = parse(Float64, s)
end

xs = 0:0.01:10
sinecurve = @lift(sin.($frequency .* xs))

lines(f[1, 1], xs, sinecurve)

f
```
\end{examplefigure}


# Toggle

A toggle with an attribute `active` that can either be true or false, to enable
or disable properties of an interactive plot.

\begin{examplefigure}{}
```julia
using GLMakie
GLMakie.activate!() # hide
fig = Figure()

ax = Axis(fig[1, 1])

toggles = [Toggle(fig, active = active) for active in [true, false]]
labels = [Label(fig, lift(x -> x ? "$l visible" : "$l invisible", t.active))
    for (t, l) in zip(toggles, ["sine", "cosine"])]

fig[1, 2] = grid!(hcat(toggles, labels), tellheight = false)

line1 = lines!(0..10, sin, color = :blue, visible = false)
line2 = lines!(0..10, cos, color = :red)

connect!(line1.visible, toggles[1].active)
connect!(line2.visible, toggles[2].active)

fig
```
\end{examplefigure}# Button

\begin{examplefigure}{}
```julia
using GLMakie
GLMakie.activate!() # hide
fig = Figure()

ax = Axis(fig[1, 1])
fig[2, 1] = buttongrid = GridLayout(tellwidth = false)

counts = Observable([1, 4, 3, 7, 2])

buttonlabels = [lift(x -> "Count: $(x[i])", counts) for i in 1:5]

buttons = buttongrid[1, 1:5] = [Button(fig, label = l) for l in buttonlabels]

for i in 1:5
    on(buttons[i].clicks) do n
        counts[][i] += 1
        notify(counts)
    end
end

barplot!(counts, color = cgrad(:Spectral)[LinRange(0, 1, 5)])
ylims!(ax, 0, 20)

fig
```
\end{examplefigure}


# Label

A Label is text within a rectangular boundingbox.
The `halign` and `valign` attributes always refer to unrotated horizontal and vertical.
This is different from `text`, where alignment is relative to text flow direction.

A Label's size is known, so if `tellwidth` and `tellheight` are set to `true` (the default values) a GridLayout with `Auto` column or row sizes can shrink to fit.

\begin{examplefigure}{svg = true}
```julia
using CairoMakie
CairoMakie.activate!() # hide

fig = Figure()

fig[1:2, 1:3] = [Axis(fig) for _ in 1:6]

supertitle = Label(fig[0, :], "Six plots", textsize = 30)

sideinfo = Label(fig[2:3, 0], "This text is vertical", rotation = pi/2)

fig
```
\end{examplefigure}

Justification and lineheight of a label can be controlled just like with normal text.

\begin{examplefigure}{svg = true}
```julia
using CairoMakie
CairoMakie.activate!() # hide

f = Figure()

Label(f[1, 1],
    "Left Justified\nMultiline\nLabel\nLineheight 0.9",
    justification = :left,
    lineheight = 0.9
)
Label(f[1, 2],
    "Center Justified\nMultiline\nLabel\nLineheight 1.1",
    justification = :center,
    lineheight = 1.1
)
Label(f[1, 3],
    "Right Justified\nMultiline\nLabel\nLineheight 1.3",
    justification = :right,
    lineheight = 1.3
)

f
```
\end{examplefigure}

# Axis

## Creating an Axis

The `Axis` is a 2D axis that works well with automatic layouts.
Here's how you create one

\begin{examplefigure}{svg = true}
```julia
using CairoMakie
CairoMakie.activate!() # hide

f = Figure()

ax = Axis(f[1, 1], xlabel = "x label", ylabel = "y label",
    title = "Title")

f
```
\end{examplefigure}
## Plotting into an Axis

You can use all the normal mutating 2D plotting functions with an `Axis`.
These functions return the created plot object.
Omitting the `ax` argument plots into the `current_axis()`, which is usually the axis that was last created.

\begin{examplefigure}{svg = true}
```julia
lineobject = lines!(ax, 0..10, sin, color = :red)
scatobject = scatter!(0:0.5:10, cos, color = :orange)

f
```
\end{examplefigure}
## Deleting plots

You can delete a plot object directly via `delete!(ax, plotobj)`.
You can also remove all plots with `empty!(ax)`.

\begin{examplefigure}{svg = true}
```julia
using CairoMakie
CairoMakie.activate!() # hide
Makie.inline!(true) # hide

f = Figure()

axs = [Axis(f[1, i]) for i in 1:3]

scatters = map(axs) do ax
    [scatter!(ax, 0:0.1:10, x -> sin(x) + i) for i in 1:3]
end

delete!(axs[2], scatters[2][2])
empty!(axs[3])

f
```
\end{examplefigure}
## Setting Axis limits and reversing axes

You can set axis limits with the functions `xlims!`, `ylims!` or `limits!`. The
numbers are meant in the order left right for `xlims!`, and bottom top for `ylims!`.
Therefore, if the second number is smaller than the first, the respective axis
will reverse. You can manually reverse an axis by setting `ax.xreversed = true` or
`ax.yreversed = true`.

Note that if you enforce an aspect ratio between x-axis and y-axis using `autolimitaspect`,
the values you set with these functions will probably not be exactly what you get,
but they will be changed to fit the chosen ratio.

\begin{examplefigure}{svg = true}
```julia
using CairoMakie
CairoMakie.activate!() # hide
Makie.inline!(true) # hide

f = Figure()

axes = [Axis(f[i, j]) for j in 1:3, i in 1:2]

for (i, ax) in enumerate(axes)
    ax.title = "Axis $i"
    poly!(ax, Point2f[(9, 9), (3, 1), (1, 3)],
        color = cgrad(:inferno, 6, categorical = true)[i])
end

xlims!(axes[1], [0, 10]) # as vector
xlims!(axes[2], 10, 0) # separate, reversed
ylims!(axes[3], 0, 10) # separate
ylims!(axes[4], (10, 0)) # as tuple, reversed
limits!(axes[5], 0, 10, 0, 10) # x1, x2, y1, y2
limits!(axes[6], BBox(0, 10, 0, 10)) # as rectangle

f
```
\end{examplefigure}
### Setting half-automatic limits

You can set half limits by either giving one argument as `nothing` or by using the keyword syntax where only `low` or `high` is given.

\begin{examplefigure}{svg = true}
```julia
using CairoMakie
CairoMakie.activate!() # hide
Makie.inline!(true) # hide

f = Figure()

data = rand(100, 2) .* 0.7 .+ 0.15

Axis(f[1, 1], title = "xlims!(nothing, 1)")
scatter!(data)
xlims!(nothing, 1)

Axis(f[1, 2], title = "xlims!(low = 0)")
scatter!(data)
xlims!(low = 0)

Axis(f[2, 1], title = "ylims!(0, nothing)")
scatter!(data)
ylims!(0, nothing)

Axis(f[2, 2], title = "ylims!(high = 1)")
scatter!(data)
ylims!(high = 1)

f
```
\end{examplefigure}
This also works when specifying limits directly, such as `Axis(..., limits = (nothing, 1, 2, nothing))`.

### Auto-reset behavior

When you create a new plot in an axis, `reset_limits!(ax)` is called, which adjusts the limits to the new bounds.
If you have previously set limits with `limits!`, `xlims!` or `ylims!`, these limits are not overridden by the new plot. If you want to override the manually set limits, call `autolimits!(ax)` to compute completely new limits from the axis content.

The user-defined limits are stored in `ax.limits`. This can either be a tuple with two entries, where each entry can be either `nothing` or a tuple with numbers `(low, high)`.It can also be a tuple with four numbers `(xlow, xhigh, ylow, yhigh)`. You can pass this directly when creating a new axis. The same observable `limits` is also set using `limits!`, `xlims!` and `ylims!`, or reset to `(nothing, nothing)` using `autolimits!`.

\begin{examplefigure}{svg = true}
```julia
using CairoMakie
CairoMakie.activate!() # hide
Makie.inline!(true) # hide

f = Figure()

lines(f[1, 1], 0..10, sin)
lines(f[1, 2], 0..10, sin, axis = (limits = (0, 10, -1, 1),))

f
```
\end{examplefigure}
## Modifying ticks

To control ticks, you can set the axis attributes `xticks/yticks` and `xtickformat/ytickformat`.

You can overload one or more of these three functions to implement custom ticks:

```julia
tickvalues, ticklabels = MakieLayout.get_ticks(ticks, scale, formatter, vmin, vmax)
tickvalues = MakieLayout.get_tickvalues(ticks, vmin, vmax)
ticklabels = MakieLayout.get_ticklabels(formatter, tickvalues)
```

If you overload `get_ticks`, you have to compute both tickvalues and ticklabels directly as a vector of floats and strings, respectively.
Otherwise the result of `get_tickvalues` is passed to `get_ticklabels` by default.
The limits of the respective axis are passed as `vmin` and `vmax`.

A couple of behaviors are implemented by default.
You can specify static ticks by passing an iterable of numbers.
You can also pass a tuple with tick values and tick labels directly, bypassing the formatting step.

As a third option you can pass a function taking minimum and maximum axis value as arguments and returning either a vector of tickvalues which are then passed to the current formatter, or a tuple with tickvalues and ticklabels which are then used directly.

For formatting, you can pass a function which takes a vector of numbers and outputs a vector of strings.
You can also pass a format string which is passed to `Formatting.format` from [Formatting.jl](https://github.com/JuliaIO/Formatting.jl), where you can mix the formatted numbers with other text like in `"{:.2f}ms"`.

### Predefined ticks

The default tick type is `LinearTicks(n)`, where `n` is the target number of ticks which the algorithm tries to return.

\begin{examplefigure}{svg = true}
```julia
using CairoMakie
CairoMakie.activate!() # hide
Makie.inline!(true) # hide

fig = Figure()
for (i, n) in enumerate([2, 5, 9])
    lines(fig[i, 1], 0..20, sin, axis = (xticks = LinearTicks(n),))
end
fig
```
\end{examplefigure}
There's also `WilkinsonTicks` which uses the alternative Wilkinson algorithm.

`MultiplesTicks` can be used when an axis should be marked at multiples of a certain number.
A common scenario is plotting a trigonometric function which should be marked at pi intervals.

\begin{examplefigure}{svg = true}
```julia
using CairoMakie
CairoMakie.activate!() # hide
Makie.inline!(true) # hide

lines(0..20, sin, axis = (xticks = MultiplesTicks(4, pi, "π"),))
```
\end{examplefigure}
Here are a couple of examples that show off different settings for ticks and formats.

\begin{examplefigure}{svg = true}
```julia
using CairoMakie
CairoMakie.activate!() # hide
Makie.inline!(true) # hide

f = Figure()

axes = [Axis(f[i, j]) for i in 1:2, j in 1:2]

xs = LinRange(0, 2pi, 50)
for (i, ax) in enumerate(axes)
    ax.title = "Axis $i"
    lines!(ax, xs, sin.(xs))
end

axes[1].xticks = 0:6

axes[2].xticks = 0:pi:2pi
axes[2].xtickformat = xs -> ["$(x/pi)π" for x in xs]

axes[3].xticks = (0:pi:2pi, ["start", "middle", "end"])

axes[4].xticks = 0:pi:2pi
axes[4].xtickformat = "{:.2f}ms"
axes[4].xlabel = "Time"

f
```
\end{examplefigure}
## Minor ticks and grids

You can show minor ticks and grids by setting `x/yminorticksvisible = true` and `x/yminorgridvisible = true` which are off by default.
You can set size, color, width, align etc. like for the normal ticks, but there are no labels.
The `x/yminorticks` attributes control how minor ticks are computed given major ticks and axis limits.
For that purpose you can create your own minortick type and overload `MakieLayout.get_minor_tickvalues(minorticks, tickvalues, vmin, vmax)`.

The default minor tick type is `IntervalsBetween(n, mirror = true)` where `n` gives the number of intervals each gap between major ticks is divided into with minor ticks, and `mirror` decides if outside of the major ticks there are more minor ticks with the same intervals as the adjacent gaps.

\begin{examplefigure}{svg = true}
```julia
using CairoMakie
CairoMakie.activate!() # hide
Makie.inline!(true) # hide

theme = Attributes(
    Axis = (
        xminorticksvisible = true,
        yminorticksvisible = true,
        xminorgridvisible = true,
        yminorgridvisible = true,
    )
)

fig = with_theme(theme) do
    fig = Figure()
    axs = [Axis(fig[fldmod1(n, 2)...],
        title = "IntervalsBetween($(n+1))",
        xminorticks = IntervalsBetween(n+1),
        yminorticks = IntervalsBetween(n+1)) for n in 1:4]
    fig
end

fig
```
\end{examplefigure}

Minor ticks can also be given as an `AbstractVector` of real numbers.

\begin{examplefigure}{svg = true}
```julia
using CairoMakie
CairoMakie.activate!() # hide
Makie.inline!(true) # hide

lines(1..10, sin, axis = (
    yminorgridvisible = true,
    yminorticksvisible = true,
    yminorticks = -0.9:0.1:0.9,
    yticks = [-1, 1],
))
```
\end{examplefigure}


## Hiding Axis spines and decorations

You can hide all axis elements manually, by setting their specific visibility attributes to `false`, like
`xticklabelsvisible`, but that can be tedious. There are a couple of convenience functions for this.

To hide spines, you can use `hidespines!`.

\begin{examplefigure}{svg = true}
```julia
using CairoMakie
CairoMakie.activate!() # hide
Makie.inline!(true) # hide

f = Figure()

ax1 = Axis(f[1, 1], title = "Axis 1")
ax2 = Axis(f[1, 2], title = "Axis 2")

hidespines!(ax1)
hidespines!(ax2, :t, :r) # only top and right

f
```
\end{examplefigure}
To hide decorations, you can use `hidedecorations!`, or the specific `hidexdecorations!` and `hideydecorations!`.
When hiding, you can set `label = false`, `ticklabels = false`, `ticks = false`, `grid = false`, `minorgrid = false` or `minorticks = false` as keyword
arguments if you want to keep those elements.
It's common, e.g., to hide everything but the grid lines in facet plots.

\begin{examplefigure}{svg = true}
```julia
using CairoMakie
CairoMakie.activate!() # hide
Makie.inline!(true) # hide

f = Figure()

ax1 = Axis(f[1, 1], title = "Axis 1")
ax2 = Axis(f[1, 2], title = "Axis 2")
ax3 = Axis(f[1, 3], title = "Axis 3")

hidedecorations!(ax1)
hidexdecorations!(ax2, grid = false)
hideydecorations!(ax3, ticks = false)

f
```
\end{examplefigure}

## Trimmed spines

The attributes `xtrimspine` and `ytrimspine` can be used to limit the respective spines to the range of the outermost major ticks.

\begin{examplefigure}{svg = true}
```julia
using CairoMakie
CairoMakie.activate!() # hide
Makie.inline!(true) # hide

hist(randn(100) ./ 4 .+ 5,
    strokewidth = 1,
    strokecolor = :black,
    axis = (
        xtrimspine = true,
        ytrimspine = true,
        topspinevisible = false,
        rightspinevisible = false,
        title = "Trimmed spines",
        xgridvisible = false,
        ygridvisible = false,
    )
)
```
\end{examplefigure}

## Log scales and other axis scales

The two attributes `xscale` and `yscale`, which by default are set to `identity`, can be used to project the data in a nonlinear way, in addition to the linear zoom that the limits provide.

Take care that the axis limits always stay inside the limits appropriate for the chosen scaling function, for example, `log` functions fail for values `x <= 0`, `sqrt` for `x < 0`, etc.

\begin{examplefigure}{svg = true}
```julia
using CairoMakie
CairoMakie.activate!() # hide
Makie.inline!(true) # hide

data = LinRange(0.01, 0.99, 200)

f = Figure(resolution = (800, 800))

for (i, scale) in enumerate([identity, log10, log2, log, sqrt, Makie.logit])

    row, col = fldmod1(i, 2)
    Axis(f[row, col], yscale = scale, title = string(scale),
        yminorticksvisible = true, yminorgridvisible = true,
        yminorticks = IntervalsBetween(8))

    lines!(data, color = :blue)
end

f
```
\end{examplefigure}
### Pseudolog and symlog scales

Some plotting functions, like barplots or density plots, have offset parameters which are usually zero, which you have to set to some non-zero value explicitly so they work in `log` axes.

\begin{examplefigure}{svg = true}
```julia
using CairoMakie
CairoMakie.activate!() # hide
Makie.inline!(true) # hide

processors = ["VAX-11/780", "Sun-4/260", "PowerPC 604",
    "Alpha 21164", "Intel Pentium III", "Intel Xeon"]
relative_speeds = [1, 9, 117, 280, 1779, 6505]

barplot(relative_speeds, fillto = 0.5,
    axis = (yscale = log10, ylabel ="relative speed",
        xticks = (1:6, processors), xticklabelrotation = pi/8))

ylims!(0.5, 10000)
current_figure()
```
\end{examplefigure}
Another option are pseudolog and symlog scales.
Pseudolog is similar to log, but modified in order to work for zero and for negative values.
The `pseudolog10` function is defined as `sign(x) * log10(abs(x) + 1)`.

Another option for symmetric log scales including zero is the symmetric log scale `Symlog10`, which combines a normal log scale with a linear scale between two boundary values around zero.

\begin{examplefigure}{svg = true}
```julia
using CairoMakie
CairoMakie.activate!() # hide
Makie.inline!(true) # hide

f = Figure(resolution = (800, 700))

lines(f[1, 1], -100:0.1:100, axis = (
    yscale = Makie.pseudolog10,
    title = "Pseudolog scale",
    yticks = [-100, -10, -1, 0, 1, 10, 100]))

lines(f[2, 1], -100:0.1:100, axis = (
    yscale = Makie.Symlog10(10.0),
    title = "Symlog10 with linear scaling between -10 and 10",
    yticks = [-100, -10, 0, 10, 100]))

f
```
\end{examplefigure}

## Controlling Axis aspect ratios

If you're plotting images, you might want to force a specific aspect ratio
of an axis, so that the images are not stretched. The default is that an axis
uses all of the available space in the layout. You can use `AxisAspect` and
`DataAspect` to control the aspect ratio. For example, `AxisAspect(1)` forces a
square axis and `AxisAspect(2)` results in a rectangle with a width of two
times the height.
`DataAspect` uses the currently chosen axis limits and brings the axes into the
same aspect ratio. This is the easiest to use with images.
A different aspect ratio can only reduce the axis space that is being used, also
it necessarily has to break the layout a little bit.

\begin{examplefigure}{svg = true}
```julia
using CairoMakie
using FileIO
using Random # hide
Random.seed!(1) # hide

f = Figure()

axes = [Axis(f[i, j]) for i in 1:2, j in 1:3]
tightlimits!.(axes)

img = rotr90(load(assetpath("cow.png")))

for ax in axes
    image!(ax, img)
end

axes[1, 1].title = "Default"

axes[1, 2].title = "DataAspect"
axes[1, 2].aspect = DataAspect()

axes[1, 3].title = "AxisAspect(418/348)"
axes[1, 3].aspect = AxisAspect(418/348)

axes[2, 1].title = "AxisAspect(1)"
axes[2, 1].aspect = AxisAspect(1)

axes[2, 2].title = "AxisAspect(2)"
axes[2, 2].aspect = AxisAspect(2)

axes[2, 3].title = "AxisAspect(2/3)"
axes[2, 3].aspect = AxisAspect(2/3)

f
```
\end{examplefigure}

## Controlling data aspect ratios

If you want the content of an axis to adhere to a certain data aspect ratio, there is
another way than forcing the aspect ratio of the whole axis to be the same, and
possibly breaking the layout. This works via the axis attribute `autolimitaspect`.
It can either be set to `nothing` which means the data limits can have any arbitrary
aspect ratio. Or it can be set to a number, in which case the targeted limits of the
axis (that are computed by `autolimits!`) are enlarged to have the correct aspect ratio.

You can see the different ways to get a plot with an unstretched circle, using
different ways of setting aspect ratios, in the following example.

```julia:video
using CairoMakie
using Animations



# scene setup for animation
###########################################################

container_scene = Scene(camera = campixel!, resolution = (1200, 1200))

t = Observable(0.0)

a_width = Animation([1, 7], [1200.0, 800], sineio(n=2, yoyo=true, postwait=0.5))
a_height = Animation([2.5, 8.5], [1200.0, 800], sineio(n=2, yoyo=true, postwait=0.5))

scene_area = lift(t) do t
    Recti(0, 0, round(Int, a_width(t)), round(Int, a_height(t)))
end

scene = Scene(container_scene, scene_area, camera = campixel!)

rect = poly!(scene, scene_area, color=RGBf(0.97, 0.97, 0.97), strokecolor=:transparent, strokewidth=0)

outer_layout = GridLayout(scene, alignmode = Outside(30))


# example begins here
###########################################################

layout = outer_layout[1, 1] = GridLayout()

titles = ["aspect enforced\nvia layout", "axis aspect\nset directly", "no aspect enforced", "data aspect conforms\nto axis size"]
axs = layout[1:2, 1:2] = [Axis(scene, title = t) for t in titles]

for a in axs
    lines!(a, Circle(Point2f(0, 0), 100f0))
end

rowsize!(layout, 1, Fixed(400))
# force the layout cell [1, 1] to be square
colsize!(layout, 1, Aspect(1, 1))

axs[2].aspect = 1
axs[4].autolimitaspect = 1

rects = layout[1:2, 1:2] = [Box(scene, color = (:black, 0.05),
    strokecolor = :transparent) for _ in 1:4]

record(container_scene, "example_circle_aspect_ratios.mp4", 0:1/30:9; framerate=30) do ti
    t[] = ti
end
nothing # hide
```

\video{example_circle_aspect_ratios}


## Linking axes

You can link axes to each other. Every axis simply keeps track of a list of other
axes which it updates when it is changed itself. You can link x and y dimensions
separately.

\begin{examplefigure}{svg = true}
```julia
using CairoMakie

f = Figure()

ax1 = Axis(f[1, 1])
ax2 = Axis(f[1, 2])
ax3 = Axis(f[2, 2])

linkyaxes!(ax1, ax2)
linkxaxes!(ax2, ax3)

ax1.title = "y linked"
ax2.title = "x & y linked"
ax3.title = "x linked"

for (i, ax) in enumerate([ax1, ax2, ax3])
    lines!(ax, 1:10, 1:10, color = "green")
    if i != 1
        lines!(ax, 11:20, 1:10, color = "red")
    end
    if i != 3
        lines!(ax, 1:10, 11:20, color = "blue")
    end
end

f
```
\end{examplefigure}
## Changing x and y axis position

By default, the x axis is at the bottom, and the y axis at the left side.
You can change this with the attributes `xaxisposition = :top` and `yaxisposition = :right`.

\begin{examplefigure}{svg = true}
```julia
using CairoMakie

f = Figure()

for i in 1:2, j in 1:2
    Axis(
        f[i, j],
        limits = (0, 5, 0, 5),
        xaxisposition = (i == 1 ? :top : :bottom),
        yaxisposition = (j == 1 ? :left : :right))
end

f
```
\end{examplefigure}
## Creating a twin axis

There is currently no dedicated function to do this, but you can simply add an Axis on top of another, then hide everything but the second axis.

Here's an example how to do this with a second y axis on the right.

\begin{examplefigure}{svg = true}
```julia
using CairoMakie

f = Figure()

ax1 = Axis(f[1, 1], yticklabelcolor = :blue)
ax2 = Axis(f[1, 1], yticklabelcolor = :red, yaxisposition = :right)
hidespines!(ax2)
hidexdecorations!(ax2)

lines!(ax1, 0..10, sin, color = :blue)
lines!(ax2, 0..10, x -> 100 * cos(x), color = :red)

f
```
\end{examplefigure}
## Axis interaction

An Axis has a couple of predefined interactions enabled.

### Scroll zoom

You can zoom in an axis by scrolling in and out.
If you press x or y while scrolling, the zoom movement is restricted to that dimension.
These keys can be changed with the attributes `xzoomkey` and `yzoomkey`.
You can also restrict the zoom dimensions all the time by setting the axis attributes `xzoomlock` or `yzoomlock` to `true`.

### Drag pan

You can pan around the axis by right-clicking and dragging.
If you press x or y while panning, the pan movement is restricted to that dimension.
These keys can be changed with the attributes `xpankey` and `ypankey`.
You can also restrict the pan dimensions all the time by setting the axis attributes `xpanlock` or `ypanlock` to `true`.

### Limit reset

You can reset the limits with `ctrl + leftclick`. This is the same as doing `reset_limits!(ax)`. This sets the limits back to the values stored in `ax.limits`, and if they are `nothing`, computes them automatically. If you have previously called `limits!`, `xlims!` or `ylims!`, these settings therefore stay intact when doing a limit reset.

You can alternatively press `ctrl + shift + leftclick`, which is the same as calling `autolimits!(ax)`.
This function ignores previously set limits and computes them all anew given the axis content.

### Rectangle selection zoom

Left-click and drag zooms into the selected rectangular area.
If you press x or y while panning, only the respective dimension is affected.
You can also restrict the selection zoom dimensions all the time by setting the axis attributes `xrectzoom` or `yrectzoom` to `true`.

### Custom interactions

The interaction system is an additional abstraction upon Makie's low-level event system to make it easier to quickly create your own interaction patterns.

#### Registering and deregistering interactions

To register a new interaction, call `register_interaction!(ax, name::Symbol, interaction)`.
The `interaction` argument can be of any type.

To remove an existing interaction completely, call `deregister_interaction!(ax, name::Symbol)`.
You can check which interactions are currently active by calling `interactions(ax)`.

#### Activating and deactivating interactions

Often, you don't want to remove an interaction entirely but only disable it for a moment, then reenable it again.
You can use the functions `activate_interaction!(ax, name::Symbol)` and `deactivate_interaction!(ax, name::Symbol)` for that.

#### `Function` interaction

If `interaction` is a `Function`, it should accept two arguments, which correspond to an event and the axis.
This function will then be called whenever the axis generates an event.

Here's an example of such a function. Note that we use the special dispatch signature for Functions that allows to use the `do`-syntax:

```julia
register_interaction!(ax, :my_interaction) do event::MouseEvent, axis
    if event.type === MouseEventTypes.leftclick
        println("You clicked on the axis!")
    end
end
```

As you can see, it's possible to restrict the type parameter of the event argument.
Choices are one of `MouseEvent`, `KeysEvent` or `ScrollEvent` if you only want to handle a specific class.
Your function can also have multiple methods dealing with each type.

#### Custom object interaction

The function option is most suitable for interactions that don't involve much state.
A more verbose but flexible option is available.
For this, you define a new type which typically holds all the state variables you're interested in.

Whenever the axis generates an event, it calls `process_interaction(interaction, event, axis)` on all
stored interactions.
By defining `process_interaction` for specific types of interaction and event, you can create more complex interaction patterns.

Here's an example with simple state handling where we allow left clicks while l is pressed, and right clicks while r is pressed:

```julia
mutable struct MyInteraction
    allow_left_click::Bool
    allow_right_click::Bool
end

function MakieLayout.process_interaction(interaction::MyInteraction, event::MouseEvent, axis)
    if interaction.use_left_click && event.type === MouseEventTypes.leftclick
        println("Left click in correct mode")
    end
    if interaction.allow_right_click && event.type === MouseEventTypes.rightclick
        println("Right click in correct mode")
    end
end

function MakieLayout.process_interaction(interaction::MyInteraction, event::KeysEvent, axis)
    interaction.allow_left_click = Keyboard.l in event.keys
    interaction.allow_right_click = Keyboard.r in event.keys
end

register_interaction!(ax, :left_and_right, MyInteraction(false, false))
```

#### Setup and cleanup

Some interactions might have more complex state involving plot objects that need to be setup or removed.
For those purposes, you can overload the methods `registration_setup!(parent, interaction)` and `deregistration_cleanup!(parent, interaction)` which are called during registration and deregistration, respectively.
# Axis3

## Viewing angles

The two attributes `azimuth` and `elevation` control the angles from which the plots are viewed.

\begin{examplefigure}{}
```julia
using GLMakie
using FileIO
GLMakie.activate!() # hide
Makie.inline!(true) # hide

f = Figure()

brain = load(assetpath("brain.stl"))
colors = [tri[1][2] for tri in brain for i in 1:3]

azimuths = [0, 0.2pi, 0.4pi]
elevations = [-0.2pi, 0, 0.2pi]

for (i, elevation) in enumerate(elevations)
    for (j, azimuth) in enumerate(azimuths)
        ax = Axis3(f[i, j], aspect = :data,
        title = "elevation = $(round(elevation/pi, digits = 2))π\nazimuth = $(round(azimuth/pi, digits = 2))π",
        elevation = elevation, azimuth = azimuth,
        protrusions = (0, 0, 0, 40))

        hidedecorations!(ax)
        mesh!(brain, color = colors, colormap = :thermal)
    end
end

f
```
\end{examplefigure}
## Data aspects and view mode

The attributes `aspect` and `viewmode` both influence the apparent relative scaling of the three axes.

### `aspect`

The `aspect` changes how long each axis is relative to the other two.

If you set it to `:data`, the axes will be scaled according to their lengths in data space.
The visual result is that objects with known real-world dimensions look correct and not squished.

You can also set it to a three-tuple, where each number gives the relative length of that axis vs the others.

\begin{examplefigure}{}
```julia
using GLMakie
using FileIO
GLMakie.activate!() # hide
Makie.inline!(true) # hide

f = Figure()

brain = load(assetpath("brain.stl"))

aspects = [:data, (1, 1, 1), (1, 2, 3), (3, 2, 1)]

for (i, aspect) in enumerate(aspects)
    ax = Axis3(f[fldmod1(i, 2)...], aspect = aspect, title = "$aspect")
    mesh!(brain, color = :bisque)
end

f
```
\end{examplefigure}
### `viewmode`

The `viewmode` changes how the final projection is adjusted to fit the axis into its scene.

The default is `:fitzoom`, which scales the final projection evenly, so that the farthest corner of the axis goes right up to the scene boundary.
If you rotate an axis with this mode, the apparent size will shrink and grow depending on the viewing angles, but the plot objects will never look skewed relative to their `aspect`.

The next option `:fit` is like `:fitzoom`, but without the zoom component.
The axis is scaled so that no matter what the viewing angles are, the axis does not clip the scene boundary and its apparent size doesn't change, even though this makes less efficient use of the available space.
You can imagine a sphere around the axis, which is zoomed right up until it touches the scene boundary.

The last option is `:stretch`.
In this mode, scaling in both x and y direction is applied to fit the axis right into its scene box.
Be aware that this mode can skew the axis a lot and doesn't keep the `aspect` intact.
On the other hand, it uses the available space most efficiently.

\begin{examplefigure}{}
```julia
using GLMakie
GLMakie.activate!() # hide
Makie.inline!(true) # hide

f = Figure()

r = LinRange(-1, 1, 100)
cube = [(x.^2 + y.^2 + z.^2) for x = r, y = r, z = r]
cube_with_holes = cube .* (cube .> 1.4)

viewmodes = [:fitzoom, :fit, :stretch]

for (j, viewmode) in enumerate(viewmodes)
    for (i, azimuth) in enumerate([1.1, 1.275, 1.45] .* pi)
        ax = Axis3(f[i, j], aspect = :data,
            azimuth = azimuth,
            viewmode = viewmode, title = "$viewmode")
        hidedecorations!(ax)
        ax.protrusions = (0, 0, 0, 20)
        volume!(cube_with_holes, algorithm = :iso, isorange = 0.05, isovalue = 1.7)
    end
end

f
```
\end{examplefigure}
## Perspective or orthographic look

You can switch smoothly between an orthographic look and a perspective look using the `perspectiveness` attribute.

A value of 0 looks like an orthographic projection (it is only approximate to a real one) while 1 gives a quite strong perspective look.

\begin{examplefigure}{}
```julia
using GLMakie
GLMakie.activate!() # hide
Makie.inline!(true) # hide

f = Figure(resolution = (1200, 800), fontsize = 14)

xs = LinRange(0, 10, 100)
ys = LinRange(0, 10, 100)
zs = [cos(x) * sin(y) for x in xs, y in ys]

for (i, perspectiveness) in enumerate(LinRange(0, 1, 6))
    Axis3(f[fldmod1(i, 3)...], perspectiveness = perspectiveness,
        title = "$perspectiveness")

    surface!(xs, ys, zs)
end

f
```
\end{examplefigure}

# Box

A simple rectangle poly that is layoutable. This can be useful to make boxes for
facet plots or when a rectangular placeholder is needed.

\begin{examplefigure}{}
```julia
using CairoMakie
using ColorSchemes
CairoMakie.activate!() # hide

fig = Figure()

rects = fig[1:4, 1:6] = [
    Box(fig, color = c)
    for c in get.(Ref(ColorSchemes.rainbow), (0:23) ./ 23)]

fig
```
\end{examplefigure}


# LScene

If you need a normal Makie scene in a layout, for example for 3D plots, you have
to use `LScene` right now. It's just a wrapper around the normal `Scene` that
makes it layoutable. The underlying Scene is accessible via the `scene` field.
You can plot into the `LScene` directly, though.

You can pass keyword arguments to the underlying `Scene` object to the `scenekw` keyword.
Currently, it can be necessary to pass a couple of attributes explicitly to make sure they are not inherited from the main scene.
To see what parameters are applicable, have a look at the [scene docs](/documentation/Scenes)

\begin{examplefigure}{}
```julia
using GLMakie
GLMakie.activate!()

fig = Figure()
pl = PointLight(Point3f(0), RGBf(20, 20, 20))
al = AmbientLight(RGBf(0.2, 0.2, 0.2))
lscene = LScene(fig[1, 1], scenekw = (lights = [pl, al], backgroundcolor=:black, clear=true), show_axis=false)
# now you can plot into lscene like you're used to
p = meshscatter!(lscene, randn(300, 3), color=:gray)
fig
```
\end{examplefigure}


# Colorbar

A Colorbar needs a colormap and a tuple of low/high limits.
The colormap's axis will then span from low to high along the visual representation of the colormap.
You can set ticks in a similar way to `Axis`.

Here's how you can create Colorbars manually.

\begin{examplefigure}{}
```julia
using CairoMakie

fig = Figure()

Axis(fig[1, 1])

# vertical colorbars
Colorbar(fig[1, 2], limits = (0, 10), colormap = :viridis,
    flipaxis = false)
Colorbar(fig[1, 3], limits = (0, 5),
    colormap = cgrad(:Spectral, 5, categorical = true), size = 25)
Colorbar(fig[1, 4], limits = (-1, 1), colormap = :heat,
    highclip = :cyan, lowclip = :red, label = "Temperature")

# horizontal colorbars
Colorbar(fig[2, 1], limits = (0, 10), colormap = :viridis,
    vertical = false)
Colorbar(fig[3, 1], limits = (0, 5), size = 25,
    colormap = cgrad(:Spectral, 5, categorical = true), vertical = false)
Colorbar(fig[4, 1], limits = (-1, 1), colormap = :heat,
    label = "Temperature", vertical = false, flipaxis = false,
    highclip = :cyan, lowclip = :red)

fig
```
\end{examplefigure}
If you pass a `plotobject`, a `heatmap` or `contourf`, the Colorbar is set up automatically such that it tracks these objects' relevant attributes like `colormap`, `colorrange`, `highclip` and `lowclip`. If you want to adjust these attributes afterwards, change them in the plot object, otherwise the Colorbar and the plot object will go out of sync.

\begin{examplefigure}{}
```julia
using CairoMakie

xs = LinRange(0, 20, 50)
ys = LinRange(0, 15, 50)
zs = [cos(x) * sin(y) for x in xs, y in ys]

fig = Figure()

ax, hm = heatmap(fig[1, 1][1, 1], xs, ys, zs)
Colorbar(fig[1, 1][1, 2], hm)

ax, hm = heatmap(fig[1, 2][1, 1], xs, ys, zs, colormap = :grays,
    colorrange = (-0.75, 0.75), highclip = :red, lowclip = :blue)
Colorbar(fig[1, 2][1, 2], hm)

ax, hm = contourf(fig[2, 1][1, 1], xs, ys, zs,
    levels = -1:0.25:1, colormap = :heat)
Colorbar(fig[2, 1][1, 2], hm, ticks = -1:0.25:1)

ax, hm = contourf(fig[2, 2][1, 1], xs, ys, zs,
    colormap = :Spectral, levels = [-1, -0.5, -0.25, 0, 0.25, 0.5, 1])
Colorbar(fig[2, 2][1, 2], hm, ticks = -1:0.25:1)

fig
```
\end{examplefigure}


# GridLayout

## Setting column and row sizes correctly

There are four different types of sizes you can give rows and columns.

### Fixed

`Fixed(scene_units)` is used to set a column or row to an absolute size, independent of its content.
This only really makes sense if there is variable width content in the column or row, that can shrink or expand to meet this size. You will probably not need `Fixed` sizes very often.

\begin{examplefigure}{}
```julia
using CairoMakie

f = Figure()

Axis(f[1, 1], title = "My column has size Fixed(400)")
Axis(f[1, 2], title = "My column has size Auto()")

colsize!(f.layout, 1, Fixed(400))
# colsize!(f.layout, 1, 400) would also work

f
```
\end{examplefigure}
### Relative

`Relative(fraction)` is used to set a column or row to a size that is a certain fraction of the available width or height.
This is useful, e.g., if you want a column to span 50% of the available width, no matter what other content is there.
In this case, you would use `Relative(1/2)`. The available width is the width of the GridLayout minus the space taken by row or column gaps including protrusions.

\begin{examplefigure}{}
```julia
using CairoMakie

f = Figure()

Axis(f[1, 1], title = "My column has size Relative(2/3)")
Axis(f[1, 2], title = "My column has size Auto()")
Colorbar(f[1, 3])

colsize!(f.layout, 1, Relative(2/3))

f
```
\end{examplefigure}

### Auto

The `Auto` size is a bit more complex to understand. It has two parameters, the Boolean `trydetermine` and the number `ratio`. The default is `Auto() == Auto(true, 1)`. This is also the default row height and column width.

A column or row that is sized `Auto(true)` tries to fit its own size to its content.

When a `GridLayout` is solved, it looks at the content in each row or column and checks if it reports a fixed size.
Many objects can report their own width or height because their content has a specific size, such as `Label`.
Other objects are often used with a user-set width or height, for example `Colorbar(scene, width = 30)`.
In this case, the `GridLayout` can also see that the colorbar has a width of 30 units and fit the column width to that value.
Objects like `Axis` on the other hand are usually not set to a specific size.

Only objects that span a *single* row or column report their width or height, respectively. If *multiple* rows or columns are spanned, it's not well defined how the space that the object needs should be distributed.

If there is more than one object with a fixed width or height in an `Auto(true)` sized column, the maximum size is used.

If a column or row is sized `Auto(false)`, fixed-size objects are ignored. It can also happen of course, that there is no fixed-size object in a row or column with `Auto(true)` size. In both these cases, columns or rows determine their size by what remains after all `Fixed`, `Relative`, `Aspect` and size-inferred `Auto(true)` columns or rows have been calculated. Each undetermined `Auto` column gets a share of the remaining space that is proportional to its `ratio` parameter.

For example, let's say there are two columns left with undetermined `Auto` size when 300 units space remain.
Column 1 has ratio `1` while column 2 has ratio `2`.
The first column will get `1 / (1 + 2) * 300 == 100` units, while the second column gets `2 / (1 + 2) * 300 == 200` units.


\begin{examplefigure}{}
```julia
using CairoMakie

f = Figure()

Axis(f[1, 1], title = "My column infers my width\nof 200 units")
Axis(f[1, 2], title = "My column gets 1/3rd\nof the remaining space")
Axis(f[1, 3], title = "My column gets 2/3rds\nof the remaining space")

colsize!(f.layout, 2, Auto(1)) # equivalent to Auto(true, 1)
colsize!(f.layout, 3, Auto(2)) # equivalent to Auto(true, 2)

f
```
\end{examplefigure}

### Aspect

This size is also a little bit trickier to understand. The syntax is `Aspect(reference, ratio)`. A column with a width of `Aspect(1, 2.0)` is set to 2.0 times the height of row 1, no matter what that height is. Therefore, the grid cell at [1, 1] will always have an aspect ratio of 2 to 1. The opposite pattern applies to rows with `Aspect` size.

Aspect sized columns or rows are very useful when you want to constrain the aspect ratio of a grid cell. For example, a plot that is always supposed to be square. Enforcing the aspect *on the layout level is better* than setting `axis.aspect = AxisAspect(1)` in most cases, because it ensures an *intact layout* where all grid cell borders are aligned visually. An `Axis` with `aspect = AxisAspect(1)` on the other hand simply shrinks so it remains square, but this will break alignment with other grid content.

\begin{examplefigure}{}
```julia
using CairoMakie

f = Figure()

Axis(f[1, 1], title = "I'm square and aligned")
Box(f[1, 2], color = (:blue, 0.1), strokecolor = :transparent)
Axis(f[1, 2], aspect = AxisAspect(1),
    title = "I'm square but break the layout.\nMy actual cell is the blue rect.")
Axis(f[2, 1])
Axis(f[2, 2])

rowsize!(f.layout, 2, Relative(2/3))
colsize!(f.layout, 1, Aspect(1, 1))

f
```
\end{examplefigure}
!!! note
    Keep in mind that if you set too many constraints on row and column sizes, a GridLayout can easily be too big or too small. It's good to have variable-width elements to fill the remaining space if you use an element with fixed size or fixed aspect ratio.

## Nesting

Grids can be nested inside other grids, and so on, to arbitrary depths. The top
grid's parent should be the scene in which the layout is placed. When you place
a grid inside another grid, that grid is automatically made its parent. Grids
also are by default set to alignmode Inside which means that the content edges
are aligned to the grid's bounding box, excluding the outer protrusions. This way,
plots in nested grids are nicely aligned along their spines.

\begin{examplefigure}{}
```julia
using CairoMakie

f = Figure()

subgl_left = GridLayout()
subgl_left[1:2, 1:2] = [Axis(f) for i in 1:2, j in 1:2]

subgl_right = GridLayout()
subgl_right[1:3, 1] = [Axis(f) for i in 1:3]

f.layout[1, 1] = subgl_left
f.layout[1, 2] = subgl_right

f
```
\end{examplefigure}

## Alignment

Here you can see the difference between the align modes Outside with and without
margins and the Inside alignmode. Only the standard Inside mode aligns the axis
spines of the contained axes nicely. The Outside mode is mostly useful for the
main GridLayout so that there some space between the window edges and the plots.
You can see that the normal axis looks the same as the one placed inside the
grid with Inside alignment, and they are both effectively aligned exactly the same.

\begin{examplefigure}{}
```julia
using CairoMakie

f = Figure(resolution = (800, 800))

Axis(f[1, 1], title = "No grid layout")
Axis(f[2, 1], title = "No grid layout")
Axis(f[3, 1], title = "No grid layout")

Axis(f[1, 2], title = "Inside", alignmode = Inside())
Axis(f[2, 2], title = "Outside", alignmode = Outside())
Axis(f[3, 2], title = "Outside(50)", alignmode = Outside(50))

[Box(f[i, 2], color = :transparent, strokecolor = :red) for i in 1:3]

f
```
\end{examplefigure}

## Spanned Placement

Elements in a grid layout can span multiple rows and columns. You can specify
them with the range syntax and colons for the full width or height. You can
also use end to specify the last row or column.

\begin{examplefigure}{}
```julia
using CairoMakie

f = Figure()

Axis(f[1, 1:2], title = "[1, 1:2]")
Axis(f[2:4, 1:2], title = "[2:4, 1:2]")
Axis(f[:, 3], title = "[:, 3]")
Axis(f[1:3, 4], title = "[1:3, 4]")
Axis(f[end, end], title = "[end, end]")

f
```
\end{examplefigure}
## Adding rows and columns by indexing

If you index outside of the current range of a grid layout, you do not get an
error. Instead, the layout automatically resizes to contain the new indices.
This is very useful if you want to iteratively build a layout, or add super or
side titles.

\begin{examplefigure}{}
```julia
using CairoMakie

f = Figure(resolution = (800, 800))

Axis(f[1, 1])
for i in 1:3
    Axis(f[:, end+1])
    Axis(f[end+1, :])
end

Label(f[0, :], text = "Super Title", textsize = 50)
Label(f[end+1, :], text = "Sub Title", textsize = 50)
Label(f[2:end-1, 0], text = "Left Text", textsize = 50,
    rotation = pi/2)
Label(f[2:end-1, end+1], text = "Right Text", textsize = 50,
    rotation = -pi/2)

f
```
\end{examplefigure}

## Trimming empty rows and columns

If you change a layout interactively and end up with unused rows or columns, `trim!`
will remove those for you.

Here we start with two axes:

\begin{examplefigure}{}
```julia
using CairoMakie

f = Figure()

ax1 = Axis(f[1, 1], title = "Axis 1")
ax2 = Axis(f[1, 2], title = "Axis 2")

f
```
\end{examplefigure}
Now we decide we'd like the second axis better if it was below the first one.
We move it two the new cell, and the old unused column is left blank.

\begin{examplefigure}{}
```julia
f[2, 1] = ax2

f
```
\end{examplefigure}
We can get rid of the unused space with `trim!`:

\begin{examplefigure}{}
```julia
trim!(f.layout)

f
```
\end{examplefigure}

## Tweaking space between rows and columns

You can use `rowgap!` and `colgap!` to change the spacing between rows or
columns respectively.

\begin{examplefigure}{}
```julia
using CairoMakie

fig = Figure()

axs = [Axis(fig[i, j]) for i in 1:3, j in 1:3]
axs[1, 1].title = "Group A"
axs[1, 2].title = "Group B.1"
axs[1, 3].title = "Group B.2"

hidedecorations!.(axs, grid=false)

colgap!(fig.layout, 1, Relative(0.15))

fig
```
\end{examplefigure}
All spaces can be changed at once by omitting the index of the gap to resize.

\begin{examplefigure}{}
```julia
rowgap!(fig.layout, 50)

fig
```
\end{examplefigure}

# contour

{{doc contour}}

### Examples

\begin{examplefigure}{}
```julia
using CairoMakie
CairoMakie.activate!() # hide
Makie.inline!(true) # hide

f = Figure()
Axis(f[1, 1])

xs = LinRange(0, 10, 100)
ys = LinRange(0, 15, 100)
zs = [cos(x) * sin(y) for x in xs, y in ys]

contour!(xs, ys, zs)

f
```
\end{examplefigure}

Omitting the `xs` and `ys` results in the indices of `zs` being used. We can also set arbitrary contour-levels using `levels`

\begin{examplefigure}{}
```julia
using CairoMakie
CairoMakie.activate!() # hide
Makie.inline!(true) # hide

f = Figure()
Axis(f[1, 1])

xs = LinRange(0, 10, 100)
ys = LinRange(0, 15, 100)
zs = [cos(x) * sin(y) for x in xs, y in ys]

contour!(zs,levels=-1:0.1:1)

f
```
\end{examplefigure}
# scatter

{{doc scatter}}

## Attributes

### Generic

- `visible::Bool = true` sets whether the plot will be rendered or not.
- `overdraw::Bool = false` sets whether the plot will draw over other plots. This specifically means ignoring depth checks in GL backends.
- `transparency::Bool = false` adjusts how the plot deals with transparency. In GLMakie `transparency = true` results in using Order Independent Transparency.
- `fxaa::Bool = false` adjusts whether the plot is rendered with fxaa (anti-aliasing). Note that scatter plots already include a different form of anti-aliasing when plotting non-image markers.
- `inspectable::Bool = true` sets whether this plot should be seen by `DataInspector`.
- `depth_shift::Float32 = 0f0` adjusts the depth value of a plot after all other transformations, i.e. in clip space, where `0 <= depth <= 1`. This only applies to GLMakie and WGLMakie and can be used to adjust render order (like a tunable overdraw).
- `model::Makie.Mat4f` sets a model matrix for the plot. This replaces adjustments made with `translate!`, `rotate!` and `scale!`.
- `color` sets the color of the plot. It can be given as a named color `Symbol` or a `Colors.Colorant`. Transparency can be included either directly as an alpha value in the `Colorant` or as an additional float in a tuple `(color, alpha)`. The color can also be set for each scattered marker by passing a `Vector` of colors or be used to index the `colormap` by passing a `Real` number or `Vector{<: Real}`.
- `colormap::Union{Symbol, Vector{<:Colorant}} = :viridis` sets the colormap that is sampled for numeric `color`s.
- `colorrange::Tuple{<:Real, <:Real}` sets the values representing the start and end points of `colormap`.
- `nan_color::Union{Symbol, <:Colorant} = RGBAf(0,0,0,0)` sets a replacement color for `color = NaN`.

### Other

- `cycle::Vector{Symbol} = [:color]` sets which attributes to cycle when creating multiple plots.
- `marker::Union{Symbol, Char, Matrix{<:Colorant}}` sets the scatter marker.
- `markersize::Union{<:Real, Vec2f} = 9` sets the size of the marker.
- `markerspace::Union{Type{Pixel}, Type{SceneSpace}} = Pixel` sets the space in which `markersize` is given. (I.e. `Pixel` units or `SceneSpace` (data) units)
- `strokewidth::Real = 0` sets the width of the outline around a marker.
- `strokecolor::Union{Symbol, <:Colorant} = :black` sets the color of the outline around a marker.
- `glowwidth::Real = 0` sets the size of a glow effect around the marker.
- `glowcolor::Union{Symbol, <:Colorant} = (:black, 0)` sets the color of the glow effect.
- `rotations::Union{Real, Billboard, Quaternion} = Billboard(0f0)` sets the rotation of the marker. A `Billboard` rotation is always around the depth axis.

## Examples

### Using x and y vectors

Scatters can be constructed by passing a list of x and y coordinates.

\begin{examplefigure}{name = "basic_scatter", svg = true}
```julia
using CairoMakie
CairoMakie.activate!() # hide
Makie.inline!(true) # hide

xs = range(0, 10, length = 30)
ys = 0.5 .* sin.(xs)

scatter(xs, ys)
```
\end{examplefigure}

### Using points

It is also possible to pass coordinates as a vector of points, which is preferred if the coordinates should be updated later, to avoid different lengths of x and y.

Attributes like `color` and `markersize` can be set in scalar or vector form.
If you pass a vector of numbers for `color`, the attribute `colorrange` which is by default automatically equal to the extrema of the color values, decides how colors are looked up in the `colormap`.

\begin{examplefigure}{svg = true}
```julia
using CairoMakie
CairoMakie.activate!() # hide
Makie.inline!(true) # hide

xs = range(0, 10, length = 30)
ys = 0.5 .* sin.(xs)
points = Point2f.(xs, ys)

scatter(points, color = 1:30, markersize = range(5, 30, length = 30),
    colormap = :thermal)
```
\end{examplefigure}

### Available markers

As markers, you can use almost any unicode character.
Currently, such glyphs are picked from the `Dejavu Sans` font, because it offers a wide range of symbols.
There is also a number of markers that can be referred to as a symbol, so that it's not necessary to find out the respective unicode character.

The backslash character examples have to be tab-completed in the REPL or editor so they are converted into unicode.

!!! note
    The scatter markers have the same sizes that the glyphs in Dejavu Sans have. This means that they are not matched in size or area. Currently, Makie does not have the option to use area matched markers, and sometimes manual adjustment might be necessary to achieve a good visual result.

\begin{examplefigure}{svg = true}
```julia
using CairoMakie
CairoMakie.activate!() # hide
Makie.inline!(true) # hide

markers_labels = [
    (:rect, ":rect"),
    (:star5, ":star5"),
    (:diamond, ":diamond"),
    (:hexagon, ":hexagon"),
    (:cross, ":cross"),
    (:xcross, ":xcross"),
    (:utriangle, ":utriangle"),
    (:dtriangle, ":dtriangle"),
    (:ltriangle, ":ltriangle"),
    (:rtriangle, ":rtriangle"),
    (:pentagon, ":pentagon"),
    (:star4, ":star4"),
    (:star8, ":star8"),
    (:vline, ":vline"),
    (:hline, ":hline"),
    (:x, ":x"),
    (:+, ":+"),
    (:circle, ":circle"),
    ('a', "'a'"),
    ('B', "'B'"),
    ('↑', "'\\uparrow'"),
    ('😄', "'\\:smile:'"),
    ('✈', "'\\:airplane:'"),
]

f = Figure()
ax = Axis(f[1, 1], yreversed = true,
    xautolimitmargin = (0.15, 0.15),
    yautolimitmargin = (0.15, 0.15)
)

for (i, (marker, label)) in enumerate(markers_labels)
    p = Point2f(fldmod1(i, 6)...)

    scatter!(p, marker = marker, markersize = 20, color = :black)
    text!(label, position = p, color = :gray70, offset = (0, 20),
        align = (:center, :bottom))
end

f
```
\end{examplefigure}

### Marker rotation

Markers can be rotated using the `rotations` attribute, which also allows to pass a vector.

\begin{examplefigure}{svg = true}
```julia
using CairoMakie
CairoMakie.activate!() # hide
Makie.inline!(true) # hide

points = [Point2f(x, y) for y in 1:10 for x in 1:10]
rotations = range(0, 2pi, length = length(points))

scatter(points, rotations = rotations, markersize = 20, marker = '↑')
```
\end{examplefigure}

### Vec markersize

You can scale x and y dimension of markers separately by passing a `Vec`.

\begin{examplefigure}{svg = true}
```julia
using CairoMakie
CairoMakie.activate!() # hide
Makie.inline!(true) # hide

f = Figure()
ax = Axis(f[1, 1])

scales = range(0.5, 1.5, length = 10)

for (i, sx) in enumerate(scales)
    for (j, sy) in enumerate(scales)
        scatter!(ax, Point2f(i, j),
            marker = '✈',
            markersize = 30 .* Vec2f(sx, sy),
            color = :black)
    end
end

f
```
\end{examplefigure}

### Marker space

By default, marker sizes do not scale relative to the data limits.
You can enable this by setting `markerspace = SceneSpace`.

\begin{examplefigure}{svg = true}
```julia
using CairoMakie
CairoMakie.activate!() # hide
Makie.inline!(true) # hide

f = Figure()
ax = Axis(f[1, 1])
limits!(ax, -10, 10, -10, 10)

scatter!(ax, Point2f(0, 0), markersize = 20, markerspace = SceneSpace,
    marker = '✈', label = "markerspace = SceneSpace")
scatter!(ax, Point2f(0, 0), markersize = 20, markerspace = Pixel,
    marker = '✈', label = "markerspace = Pixel")

axislegend(ax)

f
```
\end{examplefigure}

### Airport locations example

\begin{examplefigure}{}
```julia
using CairoMakie
using DelimitedFiles
CairoMakie.activate!() # hide
Makie.inline!(true) # hide

a = readdlm(assetpath("airportlocations.csv"))

scatter(a[1:50:end, :], marker = '✈',
    markersize = 20, color = :black)
```
\end{examplefigure}
# meshscatter

{{doc meshscatter}}

## Attributes

### Generic

- `visible::Bool = true` sets whether the plot will be rendered or not.
- `overdraw::Bool = false` sets whether the plot will draw over other plots. This specifically means ignoring depth checks in GL backends.
- `transparency::Bool = false` adjusts how the plot deals with transparency. In GLMakie `transparency = true` results in using Order Independent Transparency.
- `fxaa::Bool = true` adjusts whether the plot is rendered with fxaa (anti-aliasing).
- `inspectable::Bool = true` sets whether this plot should be seen by `DataInspector`.
- `depth_shift::Float32 = 0f0` adjusts the depth value of a plot after all other transformations, i.e. in clip space, where `0 <= depth <= 1`. This only applies to GLMakie and WGLMakie and can be used to adjust render order (like a tunable overdraw). 
- `model::Makie.Mat4f` sets a model matrix for the plot. This replaces adjustments made with `translate!`, `rotate!` and `scale!`.
- `color` sets the color of the plot. It can be given as a named color `Symbol` or a `Colors.Colorant`. Transparency can be included either directly as an alpha value in the `Colorant` or as an additional float in a tuple `(color, alpha)`. The color can also be set for each scattered mesh by passing a `Vector` of colors or be used to index the `colormap` by passing a `Real` number or `Vector{<: Real}`. 
- `colormap::Union{Symbol, Vector{<:Colorant}} = :viridis` sets the colormap that is sampled for numeric `color`s.
- `colorrange::Tuple{<:Real, <:Real}` sets the values representing the start and end points of `colormap`.
- `nan_color::Union{Symbol, <:Colorant} = RGBAf(0,0,0,0)` sets a replacement color for `color = NaN`.

### Generic 3D

- `shading = true` enables lighting.
- `diffuse::Vec3f = Vec3f(0.4)` sets how strongly the red, green and blue channel react to diffuse (scattered) light. 
- `specular::Vec3f = Vec3f(0.2)` sets how strongly the object reflects light in the red, green and blue channels.
- `shininess::Real = 32.0` sets how sharp the reflection is.
- `ssao::Bool = false` adjusts whether the plot is rendered with ssao (screen space ambient occlusion). Note that this only makes sense in 3D plots and is only applicable with `fxaa = true`.

### Other

- `cycle::Vector{Symbol} = [:color]` sets which attributes to cycle when creating multiple plots.
- `marker::Union{Symbol, GeometryBasics.GeometryPrimitive, GeometryBasics.Mesh}` sets the scattered mesh.
- `markersize::Union{<:Real, Vec3f} = 0.1` sets the scale of the mesh. This can be given as a Vector to apply to each scattered mesh individually.
- `rotations::Union{Real, Vec3f, Quaternion} = 0` sets the rotation of the mesh. A numeric rotation is around the z-axis, a `Vec3f` causes the mesh to rotate such that the the z-axis is now that vector, and a quaternion describes a general rotation. This can be given as a Vector to apply to each scattered mesh individually.


## Examples

\begin{examplefigure}{}
```julia
using GLMakie
GLMakie.activate!() # hide
Makie.inline!(true) # hide

xs = cos.(1:0.5:20)
ys = sin.(1:0.5:20)
zs = LinRange(0, 3, length(xs))

meshscatter(xs, ys, zs, markersize = 0.1, color = zs)
```
\end{examplefigure}
# mesh

{{doc mesh}}

## Attributes

### Generic

- `visible::Bool = true` sets whether the plot will be rendered or not.
- `overdraw::Bool = false` sets whether the plot will draw over other plots. This specifically means ignoring depth checks in GL backends.
- `transparency::Bool = false` adjusts how the plot deals with transparency. In GLMakie `transparency = true` results in using Order Independent Transparency.
- `fxaa::Bool = true` adjusts whether the plot is rendered with fxaa (anti-aliasing).
- `inspectable::Bool = true` sets whether this plot should be seen by `DataInspector`.
- `depth_shift::Float32 = 0f0` adjusts the depth value of a plot after all other transformations, i.e. in clip space, where `0 <= depth <= 1`. This only applies to GLMakie and WGLMakie and can be used to adjust render order (like a tunable overdraw). 
- `model::Makie.Mat4f` sets a model matrix for the plot. This replaces adjustments made with `translate!`, `rotate!` and `scale!`.
- `color` sets the color of the plot. It can be given as a named color `Symbol` or a `Colors.Colorant`. Transparency can be included either directly as an alpha value in the `Colorant` or as an additional float in a tuple `(color, alpha)`. A `Vector` of any of these can be passed to define the color per vertex. (It may be helpful to check `GeometryBasics.coordinates(my_mesh)` for this.) A `Vector{<: Real}` can also be passed to sample a colormap for each vertex. And finally, if the mesh includes uv coordinates you can pass a `Matrix` of colors to be used as a texture.
- `colormap::Union{Symbol, Vector{<:Colorant}} = :viridis` sets the colormap that is sampled for numeric `color`s.
- `colorrange::Tuple{<:Real, <:Real}` sets the values representing the start and end points of `colormap`.
- `nan_color::Union{Symbol, <:Colorant} = RGBAf(0,0,0,0)` sets a replacement color for `color = NaN`.

### Generic 3D

- `shading = true` enables lighting.
- `diffuse::Vec3f = Vec3f(0.4)` sets how strongly the red, green and blue channel react to diffuse (scattered) light. 
- `specular::Vec3f = Vec3f(0.2)` sets how strongly the object reflects light in the red, green and blue channels.
- `shininess::Real = 32.0` sets how sharp the reflection is.
- `ssao::Bool = false` adjusts whether the plot is rendered with ssao (screen space ambient occlusion). Note that this only makes sense in 3D plots and is only applicable with `fxaa = true`.


## Examples

\begin{examplefigure}{}
```julia
using GLMakie
GLMakie.activate!() # hide
Makie.inline!(true) # hide

vertices = [
    0.0 0.0;
    1.0 0.0;
    1.0 1.0;
    0.0 1.0;
]

faces = [
    1 2 3;
    3 4 1;
]

colors = [:red, :green, :blue, :orange]

scene = mesh(vertices, faces, color = colors, shading = false)
```
\end{examplefigure}

\begin{examplefigure}{}
```julia
using FileIO
using GLMakie
GLMakie.activate!() # hide
Makie.inline!(true) # hide

brain = load(assetpath("brain.stl"))

mesh(
    brain,
    color = [tri[1][2] for tri in brain for i in 1:3],
    colormap = Reverse(:Spectral),
    figure = (resolution = (1000, 1000),)
)
```
\end{examplefigure}
# image

{{doc image}}

## Attributes

### Generic

- `visible::Bool = true` sets whether the plot will be rendered or not.
- `overdraw::Bool = false` sets whether the plot will draw over other plots. This specifically means ignoring depth checks in GL backends.
- `transparency::Bool = false` adjusts how the plot deals with transparency. In GLMakie `transparency = true` results in using Order Independent Transparency.
- `fxaa::Bool = false` adjusts whether the plot is rendered with fxaa (anti-aliasing).
- `inspectable::Bool = true` sets whether this plot should be seen by `DataInspector`.
- `depth_shift::Float32 = 0f0` adjusts the depth value of a plot after all other transformations, i.e. in clip space, where `0 <= depth <= 1`. This only applies to GLMakie and WGLMakie and can be used to adjust render order (like a tunable overdraw). 
- `model::Makie.Mat4f` sets a model matrix for the plot. This replaces adjustments made with `translate!`, `rotate!` and `scale!`.
- `color` is set by the plot.
- `colormap::Union{Symbol, Vector{<:Colorant}} = [:black, :white` sets the colormap that is sampled for numeric `color`s.
- `colorrange::Tuple{<:Real, <:Real}` sets the values representing the start and end points of `colormap`.
- `nan_color::Union{Symbol, <:Colorant} = RGBAf(0,0,0,0)` sets a replacement color for `color = NaN`.

### Other

- `lowclip::Union{Nothing, Symbol, <:Colorant} = nothing` sets a color for any value below the colorrange.
- `highclip::Union{Nothing, Symbol, <:Colorant} = nothing` sets a color for any value above the colorrange.
- `interpolate::Bool = true` sets whether colors should be interpolated.


## Examples

\begin{examplefigure}{}
```julia
using CairoMakie
CairoMakie.activate!() # hide
Makie.inline!(true) # hide
using FileIO

img = load(assetpath("cow.png"))

f = Figure()

image(f[1, 1], img,
    axis = (title = "Default",))

image(f[1, 2], img,
    axis = (aspect = DataAspect(), title = "DataAspect()",))

image(f[2, 1], rotr90(img),
    axis = (aspect = DataAspect(), title = "rotr90",))

image(f[2, 2], img',
    axis = (aspect = DataAspect(), yreversed = true,
        title = "img' and reverse y-axis",))

f
```
\end{examplefigure}
# ecdfplot

{{doc ecdfplot}}

### Examples

\begin{examplefigure}{}
```julia
using CairoMakie
CairoMakie.activate!() # hide
Makie.inline!(true) # hide

f = Figure()
Axis(f[1, 1])

ecdfplot!(randn(200))

f
```
\end{examplefigure}

\begin{examplefigure}{}
```julia
using CairoMakie
CairoMakie.activate!() # hide
Makie.inline!(true) # hide

f = Figure()
Axis(f[1, 1])

x = randn(200)
ecdfplot!(x, color = (:blue, 0.3))
ecdfplot!(x, color = :red, npoints=10)

f
```
\end{examplefigure}

\begin{examplefigure}{}
```julia
using CairoMakie
CairoMakie.activate!() # hide
Makie.inline!(true) # hide

f = Figure()
Axis(f[1, 1])

x = rand(200)
w = @. x^2 * (1 - x)^2 
ecdfplot!(x)
ecdfplot!(x; weights = w, color=:orange)

f
```
\end{examplefigure}
# band

{{doc band}}

### Examples

\begin{examplefigure}{}
```julia
using CairoMakie
CairoMakie.activate!() # hide
Makie.inline!(true) # hide

f = Figure()
Axis(f[1, 1])

xs = 1:0.2:10
ys_low = -0.2 .* sin.(xs) .- 0.25
ys_high = 0.2 .* sin.(xs) .+ 0.25

band!(xs, ys_low, ys_high)
band!(xs, ys_low .- 1, ys_high .-1, color = :red)

f
```
\end{examplefigure}

\begin{examplefigure}{}
```julia
using Statistics
using CairoMakie
CairoMakie.activate!() # hide
Makie.inline!(true) # hide

f = Figure()
Axis(f[1, 1])

n, m = 100, 101
t = range(0, 1, length=m)
X = cumsum(randn(n, m), dims = 2)
X = X .- X[:, 1]
μ = vec(mean(X, dims=1)) # mean
lines!(t, μ)              # plot mean line
σ = vec(std(X, dims=1))  # stddev
band!(t, μ + σ, μ - σ)   # plot stddev band
f
```
\end{examplefigure}

\begin{examplefigure}{}
```julia
using GLMakie
GLMakie.activate!() # hide
Makie.inline!(true) # hide
lower = fill(Point3f(0,0,0), 100)
upper = [Point3f(sin(x), cos(x), 1.0) for x in range(0,2pi, length=100)]
col = repeat([1:50;50:-1:1],outer=2)
band(lower, upper, color=col, axis=(type=Axis3,))
```
\end{examplefigure}
# volume

{{doc volume}}

## Attributes

### Generic

- `visible::Bool = true` sets whether the plot will be rendered or not.
- `overdraw::Bool = false` sets whether the plot will draw over other plots. This specifically means ignoring depth checks in GL backends.
- `transparency::Bool = false` adjusts how the plot deals with transparency. In GLMakie `transparency = true` results in using Order Independent Transparency.
- `fxaa::Bool = true` adjusts whether the plot is rendered with fxaa (anti-aliasing).
- `inspectable::Bool = true` sets whether this plot should be seen by `DataInspector`.
- `depth_shift::Float32 = 0f0` adjusts the depth value of a plot after all other transformations, i.e. in clip space, where `0 <= depth <= 1`. This only applies to GLMakie and WGLMakie and can be used to adjust render order (like a tunable overdraw). 
- `model::Makie.Mat4f` sets a model matrix for the plot. This replaces adjustments made with `translate!`, `rotate!` and `scale!`.
- `colormap::Union{Symbol, Vector{<:Colorant}} = :viridis` sets the colormap that is sampled for numeric `color`s.
- `colorrange::Tuple{<:Real, <:Real}` sets the values representing the start and end points of `colormap`.
- `nan_color::Union{Symbol, <:Colorant} = RGBAf(0,0,0,0)` sets a replacement color for `color = NaN`.

### Generic 3D

- `shading = true` enables lighting.
- `diffuse::Vec3f = Vec3f(0.4)` sets how strongly the red, green and blue channel react to diffuse (scattered) light. 
- `specular::Vec3f = Vec3f(0.2)` sets how strongly the object reflects light in the red, green and blue channels.
- `shininess::Real = 32.0` sets how sharp the reflection is.
- `ssao::Bool = false` adjusts whether the plot is rendered with ssao (screen space ambient occlusion). Note that this only makes sense in 3D plots and is only applicable with `fxaa = true`.

### Other

- `algorithm::Union{Symbol, RaymarchAlgorithm} = :mip` sets the volume algorithm that is used.
- `isorange::Real = 0.05` sets the range of values picked up by the IsoValue algorithm.
- `isovalue = 0.5` sets the target value for the IsoValue algorithm.


## Examples

\begin{examplefigure}{}
```julia
using GLMakie
GLMakie.activate!() # hide
r = LinRange(-1, 1, 100)
cube = [(x.^2 + y.^2 + z.^2) for x = r, y = r, z = r]
contour(cube, alpha=0.5)
```
\end{examplefigure}

\begin{examplefigure}{}
```julia
cube_with_holes = cube .* (cube .> 1.4)
volume(cube_with_holes, algorithm = :iso, isorange = 0.05, isovalue = 1.7)
```
\end{examplefigure}

\begin{examplefigure}{}
```julia
using NIfTI
brain = niread(Makie.assetpath("brain.nii.gz")).raw
mini, maxi = extrema(brain)
normed = Float32.((brain .- mini) ./ (maxi - mini))

fig = Figure(resolution=(1000, 450))
# Make a colormap, with the first value being transparent
colormap = RGBAf.(to_colormap(:plasma), 1.0)
colormap[1] = RGBAf(0,0,0,0)
volume(fig[1, 1], normed, algorithm = :absorption, absorption=4f0, colormap=colormap, axis=(type=Axis3, title = "Absorption"))
volume(fig[1, 2], normed, algorithm = :mip, colormap=colormap, axis=(type=Axis3, title="Maximum Intensity Projection"))
fig
```
\end{examplefigure}
# qqplot and qqnorm

{{doc qqplot}}
{{doc qqnorm}}

### Examples

Test if `xs` and `ys` follow the same distribution.

\begin{examplefigure}{}
```julia
using CairoMakie
CairoMakie.activate!() # hide
Makie.inline!(true) # hide

xs = randn(100)
ys = randn(100)

qqplot(xs, ys, qqline = :identity)
```
\end{examplefigure}

Test if `ys` is normally distributed.

\begin{examplefigure}{}
```julia
using CairoMakie
CairoMakie.activate!() # hide
Makie.inline!(true) # hide

ys = 2 .* randn(100) .+ 3

qqnorm(ys, qqline = :fitrobust)
```
\end{examplefigure}
# hlines! and vlines!

{{doc hlines!}}
{{doc vlines!}}

These functions are not plot types / recipes and only work with `Axis`.

\begin{examplefigure}{}
```julia
using CairoMakie
CairoMakie.activate!() # hide

f = Figure()

ax1 = Axis(f[1, 1], title = "vlines")

lines!(ax1, 0..4pi, sin)
vlines!(ax1, [pi, 2pi, 3pi], color = :red)

ax2 = Axis(f[1, 2], title = "hlines")
hlines!(ax2, [1, 2, 3, 4], xmax = [0.25, 0.5, 0.75, 1], color = :blue)

f
```
\end{examplefigure}# text

{{doc text}}

## Attributes

### Generic

- `visible::Bool = true` sets whether the plot will be rendered or not.
- `overdraw::Bool = false` sets whether the plot will draw over other plots. This specifically means ignoring depth checks in GL backends.
- `transparency::Bool = false` adjusts how the plot deals with transparency. In GLMakie `transparency = true` results in using Order Independent Transparency.
- `fxaa::Bool = false` adjusts whether the plot is rendered with fxaa (anti-aliasing). Note that text plots already include a different form of anti-aliasing.
- `inspectable::Bool = true` sets whether this plot should be seen by `DataInspector`.
- `depth_shift::Float32 = 0f0` adjusts the depth value of a plot after all other transformations, i.e. in clip space, where `0 <= depth <= 1`. This only applies to GLMakie and WGLMakie and can be used to adjust render order (like a tunable overdraw). 
- `model::Makie.Mat4f` sets a model matrix for the plot. This replaces adjustments made with `translate!`, `rotate!` and `scale!`.
- `color` sets the color of the plot. It can be given as a named color `Symbol` or a `Colors.Colorant`. Transparency can be included either directly as an alpha value in the `Colorant` or as an additional float in a tuple `(color, alpha)`. The color can also be set for each character by passing a `Vector` of colors. 

### Other

- `align::Tuple{Union{Symbol, Real}, Union{Symbol, Real}} = (:left, :bottom)` sets the alignment of the string w.r.t. `position`. Uses `:left, :center, :right, :top, :bottom, :baseline` or fractions.
- `font::Union{String, Vector{String}} = "Dejavu Sans"` sets the font for the string or each character.
- `justification::Union{Real, Symbol} = automatic` sets the alignment of text w.r.t its bounding box. Can be `:left, :center, :right` or a fraction. Will default to the horizontal alignment in `align`.
- `position::Union{Point2f, Point3f} = Point2f(0)` sets an anchor position for text. Can also be a `Vector` of positions.
- `rotation::Union{Real, Quaternion}` rotates text around the given position.
- `textsize::Union{Real, Vec2f}` sets the size of each character.
- `space::Symbol = :screen` sets the space in which `textsize` acts. Can be `:screen` (pixel space) or `:data` (world space).
- `strokewidth::Real = 0` sets the width of the outline around a marker.
- `strokecolor::Union{Symbol, <:Colorant} = :black` sets the color of the outline around a marker.
- `glowwidth::Real = 0` sets the size of a glow effect around the marker.
- `glowcolor::Union{Symbol, <:Colorant} = (:black, 0)` sets the color of the glow effect.


## Screen space text

By default, text is drawn in screen space (`space = :screen`).
The text anchor is given in data coordinates, but the size of the glyphs is independent of data scaling.
The boundingbox of the text will include every data point or every text anchor point.
This also means that `autolimits!` might cut off your text, because the glyphs don't have a meaningful size in data coordinates (the size is independent of zoom level), and you have to take some care to manually place the text or set data limits such that it is fully visible.

You can either plot one string with one position, or a vector of strings with a vector of positions.

\begin{examplefigure}{}
```julia
using CairoMakie
CairoMakie.activate!() # hide
Makie.inline!(true) # hide

f = Figure()

Axis(f[1, 1], aspect = DataAspect(), backgroundcolor = :gray50)

scatter!(Point2f(0, 0))
text!("center", position = (0, 0), align = (:center, :center))

circlepoints = [(cos(a), sin(a)) for a in LinRange(0, 2pi, 16)[1:end-1]]
scatter!(circlepoints)
text!(
    "this is point " .* string.(1:15),
    position = circlepoints,
    rotation = LinRange(0, 2pi, 16)[1:end-1],
    align = (:right, :baseline),
    color = cgrad(:Spectral)[LinRange(0, 1, 15)]
)

f
```
\end{examplefigure}

## Data space text

For text whose dimensions are meaningful in data space, set `space = :data`.
This means that the boundingbox of the text in data coordinates will include every glyph.

\begin{examplefigure}{}
```julia
using CairoMakie
CairoMakie.activate!() # hide
Makie.inline!(true) # hide

f = Figure()
LScene(f[1, 1])

text!(
    fill("Makie", 7),
    rotation = [i / 7 * 1.5pi for i in 1:7],
    position = [Point3f(0, 0, i/2) for i in 1:7],
    color = [cgrad(:viridis)[x] for x in LinRange(0, 1, 7)],
    align = (:left, :baseline),
    textsize = 1,
    space = :data
)

f
```
\end{examplefigure}

## Justification

By default, justification of multiline text follows alignment.
Text that is left aligned is also left justified.
You can override this with the `justification` attribute.

\begin{examplefigure}{}
```julia
using CairoMakie
CairoMakie.activate!() # hide
Makie.inline!(true) # hide

scene = Scene(camera = campixel!, resolution = (800, 800))

points = [Point(x, y) .* 200 for x in 1:3 for y in 1:3]
scatter!(scene, points, marker = :circle, markersize = 10px)

symbols = (:left, :center, :right)

for ((justification, halign), point) in zip(Iterators.product(symbols, symbols), points)

    t = text!(scene, "a\nshort\nparagraph",
        color = (:black, 0.5),
        position = point,
        align = (halign, :center),
        justification = justification)

    bb = boundingbox(t)
    wireframe!(scene, bb, color = (:red, 0.2))
end

for (p, al) in zip(points[3:3:end], (:left, :center, :right))
    text!(scene, "align :" * string(al), position = p .+ (0, 80),
        align = (:center, :baseline))
end

for (p, al) in zip(points[7:9], (:left, :center, :right))
    text!(scene, "justification\n:" * string(al), position = p .+ (80, 0),
        align = (:center, :top), rotation = pi/2)
end

scene
```
\end{examplefigure}

## Offset

The offset attribute can be used to shift text away from its position.
This is especially useful with `space = :screen`, for example to place text together with barplots.
You can specify the end of the barplots in data coordinates, and then offset the text a little bit to the left.

\begin{examplefigure}{}
```julia
using CairoMakie
CairoMakie.activate!() # hide
Makie.inline!(true) # hide

f = Figure()

horsepower = [52, 78, 80, 112, 140]
cars = ["Kia", "Mini", "Honda", "Mercedes", "Ferrari"]

ax = Axis(f[1, 1], xlabel = "horse power")
tightlimits!(ax, Left())
hideydecorations!(ax)

barplot!(horsepower, direction = :x)
text!(cars, position = Point.(horsepower, 1:5), align = (:right, :center),
    offset = (-20, 0), color = :white)

f
```
\end{examplefigure}

## MathTeX

Makie can render LaTeX strings from the LaTeXStrings.jl package using [MathTeXEngine.jl](https://github.com/Kolaru/MathTeXEngine.jl/).
For example, you can pass L-strings as labels to the legend.

\begin{examplefigure}{}
```julia
using CairoMakie
CairoMakie.activate!() # hide
Makie.inline!(true) # hide

f = Figure()
ax = Axis(f[1, 1])

lines!(0..10, x -> sin(3x) / (cos(x) + 2),
    label = L"\frac{\sin(3x)}{\cos(x) + 2}")
lines!(0..10, x -> sin(x^2) / (cos(sqrt(x)) + 2),
    label = L"\frac{\sin(x^2)}{\cos(\sqrt{x}) + 2}")

Legend(f[1, 2], ax)

f
```
\end{examplefigure}
# streamplot

{{doc streamplot}}

### Examples

\begin{examplefigure}{}
```julia
using CairoMakie
CairoMakie.activate!() # hide
Makie.inline!(true) # hide

struct FitzhughNagumo{T}
    ϵ::T
    s::T
    γ::T
    β::T
end

P = FitzhughNagumo(0.1, 0.0, 1.5, 0.8)

f(x, P::FitzhughNagumo) = Point2f(
    (x[1]-x[2]-x[1]^3+P.s)/P.ϵ,
    P.γ*x[1]-x[2] + P.β
)

f(x) = f(x, P)

streamplot(f, -1.5..1.5, -1.5..1.5, colormap = :magma)
```
\end{examplefigure}
# poly

{{doc poly}}

### Examples

\begin{examplefigure}{}
```julia
using CairoMakie
CairoMakie.activate!() # hide
using Makie.GeometryBasics
Makie.inline!(true) # hide

f = Figure()
Axis(f[1, 1])

poly!(Point2f[(0, 0), (2, 0), (3, 1), (1, 1)], color = :red, strokecolor = :black, strokewidth = 1)

f
```
\end{examplefigure}

\begin{examplefigure}{}
```julia
using CairoMakie
CairoMakie.activate!() # hide
using Makie.GeometryBasics
Makie.inline!(true) # hide

f = Figure()
Axis(f[1, 1])

# polygon with hole
p = Polygon(
    Point2f[(0, 0), (2, 0), (3, 1), (1, 1)],
    [Point2f[(0.75, 0.25), (1.75, 0.25), (2.25, 0.75), (1.25, 0.75)]]
)

poly!(p, color = :blue)

f
```
\end{examplefigure}

\begin{examplefigure}{}
```julia
using CairoMakie
CairoMakie.activate!() # hide
using Makie.GeometryBasics
Makie.inline!(true) # hide

f = Figure()
Axis(f[1, 1])

# vector of shapes
poly!(
    [Rect(i, j, 0.75, 0.5) for i in 1:5 for j in 1:3],
    color = 1:15,
    colormap = :heat
)

f
```
\end{examplefigure}

\begin{examplefigure}{}
```julia
using CairoMakie
CairoMakie.activate!() # hide
using Makie.GeometryBasics
Makie.inline!(true) # hide

f = Figure()
Axis(f[1, 1], aspect = DataAspect())

# shape decomposition
poly!(Circle(Point2f(0, 0), 15f0), color = :pink)

f
```
\end{examplefigure}

\begin{examplefigure}{}
```julia
using CairoMakie
CairoMakie.activate!() # hide
using Makie.GeometryBasics
Makie.inline!(true) # hide

f = Figure()
Axis(f[1, 1]; backgroundcolor = :gray15)

# vector of polygons
ps = [Polygon(rand(Point2f, 3) .+ Point2f(i, j))
    for i in 1:5 for j in 1:10]

poly!(ps, color = rand(RGBf, length(ps)))

f
```
\end{examplefigure}
# errorbars

{{doc errorbars}}

### Examples

\begin{examplefigure}{}
```julia
using CairoMakie
CairoMakie.activate!() # hide
Makie.inline!(true) # hide

f = Figure()
Axis(f[1, 1])

xs = 0:0.5:10
ys = 0.5 .* sin.(xs)

lowerrors = fill(0.1, length(xs))
higherrors = LinRange(0.1, 0.4, length(xs))

errorbars!(xs, ys, higherrors, color = :red) # same low and high error

# plot position scatters so low and high errors can be discriminated
scatter!(xs, ys, markersize = 3, color = :black)

f
```
\end{examplefigure}

\begin{examplefigure}{}
```julia
using CairoMakie
CairoMakie.activate!() # hide
Makie.inline!(true) # hide

f = Figure()
Axis(f[1, 1])

xs = 0:0.5:10
ys = 0.5 .* sin.(xs)

lowerrors = fill(0.1, length(xs))
higherrors = LinRange(0.1, 0.4, length(xs))

errorbars!(xs, ys, lowerrors, higherrors,
    color = range(0, 1, length = length(xs)),
    whiskerwidth = 10)

# plot position scatters so low and high errors can be discriminated
scatter!(xs, ys, markersize = 3, color = :black)

f
```
\end{examplefigure}

\begin{examplefigure}{}
```julia
using CairoMakie
CairoMakie.activate!() # hide
Makie.inline!(true) # hide

f = Figure()
Axis(f[1, 1])

xs = 0:0.5:10
ys = 0.5 .* sin.(xs)

lowerrors = fill(0.1, length(xs))
higherrors = LinRange(0.1, 0.4, length(xs))

errorbars!(xs, ys, lowerrors, higherrors, whiskerwidth = 3, direction = :x)

# plot position scatters so low and high errors can be discriminated
scatter!(xs, ys, markersize = 3, color = :black)

f
```
\end{examplefigure}
# volumeslices

{{doc volumeslices}}

### Examples

\begin{examplefigure}{}
```julia
using GLMakie
GLMakie.activate!() # hide
Makie.inline!(true) # hide

fig = Figure()
ax = LScene(fig[1, 1], scenekw=(show_axis=false,))

x = LinRange(0, π, 50)
y = LinRange(0, 2π, 100)
z = LinRange(0, 3π, 150)

lsgrid = labelslidergrid!(
  fig,
  ["yz plane - x axis", "xz plane - y axis", "xy plane - z axis"],
  [1:length(x), 1:length(y), 1:length(z)]
)
fig[2, 1] = lsgrid.layout

vol = [cos(X)*sin(Y)*sin(Z) for X ∈ x, Y ∈ y, Z ∈ z]
plt = volumeslices!(ax, x, y, z, vol)

# connect sliders to `volumeslices` update methods
sl_yz, sl_xz, sl_xy = lsgrid.sliders

on(sl_yz.value) do v; plt[:update_yz][](v) end
on(sl_xz.value) do v; plt[:update_xz][](v) end
on(sl_xy.value) do v; plt[:update_xy][](v) end

set_close_to!(sl_yz, .5length(x))
set_close_to!(sl_xz, .5length(y))
set_close_to!(sl_xy, .5length(z))

# cam3d!(ax.scene, projectiontype=Makie.Orthographic)

fig
```
\end{examplefigure}
# abline!

{{doc abline!}}

The function `abline!` draws a linear function given slope and intercept values through the given `Axis`. The line always spans across the whole axis and doesn't affect the limits.

This function is not a plot type / recipe and only works with `Axis`.

\begin{examplefigure}{}
```julia
using CairoMakie
CairoMakie.activate!() # hide

fig, ax, pl = scatter(1:4)
abline!(ax, 0, 1)
abline!(ax, 0, 1.5, color = :red, linestyle=:dash, linewidth=2)
fig
```
\end{examplefigure}# linesegments

{{doc linesegments}}

## Attributes

### Generic

- `visible::Bool = true` sets whether the plot will be rendered or not.
- `overdraw::Bool = false` sets whether the plot will draw over other plots. This specifically means ignoring depth checks in GL backends.
- `transparency::Bool = false` adjusts how the plot deals with transparency. In GLMakie `transparency = true` results in using Order Independent Transparency.
- `fxaa::Bool = false` adjusts whether the plot is rendered with fxaa (anti-aliasing). Note that line plots already use a different form of anti-aliasing.
- `inspectable::Bool = true` sets whether this plot should be seen by `DataInspector`.
- `depth_shift::Float32 = 0f0` adjusts the depth value of a plot after all other transformations, i.e. in clip space, where `0 <= depth <= 1`. This only applies to GLMakie and WGLMakie and can be used to adjust render order (like a tunable overdraw). 
- `model::Makie.Mat4f` sets a model matrix for the plot. This replaces adjustments made with `translate!`, `rotate!` and `scale!`.
- `color` sets the color of the plot. It can be given as a named color `Symbol` or a `Colors.Colorant`. Transparency can be included either directly as an alpha value in the `Colorant` or as an additional float in a tuple `(color, alpha)`. The color can also be set for each point in the line by passing a `Vector` or be used to index the `colormap` by passing a `Real` number or `Vector{<: Real}`. 
- `colormap::Union{Symbol, Vector{<:Colorant}} = :viridis` sets the colormap that is sampled for numeric `color`s.
- `colorrange::Tuple{<:Real, <:Real}` sets the values representing the start and end points of `colormap`.
- `nan_color::Union{Symbol, <:Colorant} = RGBAf(0,0,0,0)` sets a replacement color for `color = NaN`.

### Other

- `cycle::Vector{Symbol} = [:color]` sets which attributes to cycle when creating multiple plots.
- `linestyle::Union{Nothing, Symbol, Vector} = nothing` sets the pattern of the line (e.g. `:solid`, `:dot`, `:dashdot`)
- `linewidth::Real = 1.5` sets the width of the line in pixel units.


## Examples

\begin{examplefigure}{}
```julia
using CairoMakie
CairoMakie.activate!() # hide
Makie.inline!(true) # hide

f = Figure()
Axis(f[1, 1])

xs = 1:0.2:10
ys = sin.(xs)

linesegments!(xs, ys)
linesegments!(xs, ys .- 1, linewidth = 5)
linesegments!(xs, ys .- 2, linewidth = 5, color = LinRange(1, 5, length(xs)))

f
```
\end{examplefigure}
# density

{{doc density}}

### Examples

\begin{examplefigure}{}
```julia
using CairoMakie
CairoMakie.activate!() # hide
Makie.inline!(true) # hide

f = Figure()
Axis(f[1, 1])

density!(randn(200))

f
```
\end{examplefigure}

\begin{examplefigure}{}
```julia
using CairoMakie
CairoMakie.activate!() # hide
Makie.inline!(true) # hide

f = Figure()
Axis(f[1, 1])

density!(randn(200), direction = :y, npoints = 10)

f
```
\end{examplefigure}

\begin{examplefigure}{}
```julia
using CairoMakie
CairoMakie.activate!() # hide
Makie.inline!(true) # hide

f = Figure()
Axis(f[1, 1])

density!(randn(200), color = (:red, 0.3),
    strokecolor = :red, strokewidth = 3, strokearound = true)

f
```
\end{examplefigure}

\begin{examplefigure}{}
```julia
using CairoMakie
CairoMakie.activate!() # hide
Makie.inline!(true) # hide

f = Figure()
Axis(f[1, 1])

vectors = [randn(1000) .+ i/2 for i in 0:5]

for (i, vector) in enumerate(vectors)
    density!(vector, offset = -i/4, color = (:slategray, 0.4),
        bandwidth = 0.1)
end

f
```
\end{examplefigure}

#### Gradients

You can color density plots with gradients by choosing `color = :x` or `:y`, depending on the `direction` attribute.

\begin{examplefigure}{}
```julia
using CairoMakie
CairoMakie.activate!() # hide
Makie.inline!(true) # hide

months = ["January", "February", "March", "April",
    "May", "June", "July", "August", "September",
    "October", "November", "December"]

f = Figure()
Axis(f[1, 1], title = "Fictive temperatures",
    yticks = ((1:12) ./ 4,  reverse(months)))

for i in 12:-1:1
    d = density!(randn(200) .- 2sin((i+3)/6*pi), offset = i / 4,
        color = :x, colormap = :thermal, colorrange = (-5, 5),
        strokewidth = 1, strokecolor = :black)
    # this helps with layering in GLMakie
    translate!(d, 0, 0, -0.1i)
end
f
```
\end{examplefigure}

Due to technical limitations, if you color the `:vertical` dimension (or :horizontal with direction = :y), only a colormap made with just two colors can currently work:

\begin{examplefigure}{}
```julia
using CairoMakie
CairoMakie.activate!() # hide
Makie.inline!(true) # hide

f = Figure()
Axis(f[1, 1])
for x in 1:5
    d = density!(x * randn(200) .+ 3x,
        color = :y, colormap = [:darkblue, :gray95])
end
f
```
\end{examplefigure}
# crossbar

{{doc crossbar}}

### Examples

\begin{examplefigure}{}
```julia
using CairoMakie
CairoMakie.activate!() # hide
Makie.inline!(true) # hide

xs = [1, 1, 2, 2, 3, 3]
ys = rand(6)
ymins = ys .- 1
ymaxs = ys .+ 1
dodge = [1, 2, 1, 2, 1, 2]

crossbar(xs, ys, ymins, ymaxs, dodge = dodge, show_notch = true)
```
\end{examplefigure}
# contourf

{{doc contourf}}

\begin{examplefigure}{}
```julia
using CairoMakie
CairoMakie.activate!() # hide
Makie.inline!(true) # hide

xs = LinRange(0, 10, 100)
ys = LinRange(0, 10, 100)
zs = [cos(x) * sin(y) for x in xs, y in ys]

f = Figure()
Axis(f[1, 1])

co = contourf!(xs, ys, zs, levels = 10)

Colorbar(f[1, 2], co)

f
```
\end{examplefigure}

\begin{examplefigure}{}
```julia
using CairoMakie
CairoMakie.activate!() # hide
Makie.inline!(true) # hide

xs = LinRange(0, 10, 100)
ys = LinRange(0, 10, 100)
zs = [cos(x) * sin(y) for x in xs, y in ys]

f = Figure()
Axis(f[1, 1])

co = contourf!(xs, ys, zs, levels = -0.75:0.25:0.5,
    extendlow = :cyan, extendhigh = :magenta)

Colorbar(f[1, 2], co)

f
```
\end{examplefigure}

\begin{examplefigure}{}
```julia
using CairoMakie
CairoMakie.activate!() # hide
Makie.inline!(true) # hide

xs = LinRange(0, 10, 100)
ys = LinRange(0, 10, 100)
zs = [cos(x) * sin(y) for x in xs, y in ys]

f = Figure()
Axis(f[1, 1])

co = contourf!(xs, ys, zs,
    levels = -0.75:0.25:0.5,
    extendlow = :auto, extendhigh = :auto)

Colorbar(f[1, 2], co)

f
```
\end{examplefigure}

#### Relative mode

Sometimes it's beneficial to drop one part of the range of values, usually towards the outer boundary.
Rather than specifying the levels to include manually, you can set the `mode` attribute
to `:relative` and specify the levels from 0 to 1, relative to the current minimum and maximum value.

\begin{examplefigure}{}
```julia
using CairoMakie
CairoMakie.activate!() # hide
Makie.inline!(true) # hide


using Makie.KernelDensity

k = kde([randn() + rand([0, 5]) for i in 1:10000, j in 1:2])

f = Figure(resolution = (800, 500))

Axis(f[1, 1], title = "Relative mode, drop lowest 10%")
contourf!(k, levels = 0.1:0.1:1, mode = :relative)

Axis(f[1, 2], title = "Normal mode")
contourf!(k, levels = 10)

f
```
\end{examplefigure}
# rangebars

{{doc rangebars}}

### Examples

\begin{examplefigure}{}
```julia
using CairoMakie
CairoMakie.activate!() # hide
Makie.inline!(true) # hide

f = Figure()
Axis(f[1, 1])

vals = -1:0.1:1
lows = zeros(length(vals))
highs = LinRange(0.1, 0.4, length(vals))

rangebars!(vals, lows, highs, color = :red)

f
```
\end{examplefigure}

\begin{examplefigure}{}
```julia
using CairoMakie
CairoMakie.activate!() # hide
Makie.inline!(true) # hide

f = Figure()
Axis(f[1, 1])

vals = -1:0.1:1
lows = zeros(length(vals))
highs = LinRange(0.1, 0.4, length(vals))

rangebars!(vals, lows, highs, color = LinRange(0, 1, length(vals)),
    whiskerwidth = 10, direction = :x)

f
```
\end{examplefigure}
# lines

{{doc lines}}

## Attributes

### Generic

- `visible::Bool = true` sets whether the plot will be rendered or not.
- `overdraw::Bool = false` sets whether the plot will draw over other plots. This specifically means ignoring depth checks in GL backends.
- `transparency::Bool = false` adjusts how the plot deals with transparency. In GLMakie `transparency = true` results in using Order Independent Transparency.
- `fxaa::Bool = false` adjusts whether the plot is rendered with fxaa (anti-aliasing). Note that line plots already use a different form of anti-aliasing.
- `inspectable::Bool = true` sets whether this plot should be seen by `DataInspector`.
- `depth_shift::Float32 = 0f0` adjusts the depth value of a plot after all other transformations, i.e. in clip space, where `0 <= depth <= 1`. This only applies to GLMakie and WGLMakie and can be used to adjust render order (like a tunable overdraw). 
- `model::Makie.Mat4f` sets a model matrix for the plot. This replaces adjustments made with `translate!`, `rotate!` and `scale!`.
- `color` sets the color of the plot. It can be given as a named color `Symbol` or a `Colors.Colorant`. Transparency can be included either directly as an alpha value in the `Colorant` or as an additional float in a tuple `(color, alpha)`. The color can also be set for each point in the line by passing a `Vector` of colors or be used to index the `colormap` by passing a `Real` number or `Vector{<: Real}`.
- `colormap::Union{Symbol, Vector{<:Colorant}} = :viridis` sets the colormap that is sampled for numeric `color`s.
- `colorrange::Tuple{<:Real, <:Real}` sets the values representing the start and end points of `colormap`.
- `nan_color::Union{Symbol, <:Colorant} = RGBAf(0,0,0,0)` sets a replacement color for `color = NaN`.

### Other

- `cycle::Vector{Symbol} = [:color]` sets which attributes to cycle when creating multiple plots.
- `linestyle::Union{Nothing, Symbol, Vector} = nothing` sets the pattern of the line (e.g. `:solid`, `:dot`, `:dashdot`)
- `linewidth::Real = 1.5` sets the width of the line in pixel units.


## Examples

\begin{examplefigure}{svg = true}
```julia
using CairoMakie
CairoMakie.activate!() # hide
Makie.inline!(true) # hide

f = Figure()
Axis(f[1, 1])

xs = 0:0.01:10
ys = 0.5 .* sin.(xs)

lines!(xs, ys)
lines!(xs, ys .- 1, linewidth = 5)
lines!(xs, ys .- 2, linewidth = 5, color = ys)
lines!(xs, ys .- 3, linestyle = :dash)

f
```
\end{examplefigure}

### Linestyles

\begin{examplefigure}{svg = true}
```julia
using CairoMakie
CairoMakie.activate!() # hide
Makie.inline!(true) # hide

f = Figure()
Axis(f[1, 1])

xs = 0:0.01:10
ys = 0.5 .* sin.(xs)

for (i, lw) in enumerate([1, 2, 3])
    lines!(xs, ys .- i/6, linestyle = nothing, linewidth = lw)
    lines!(xs, ys .- i/6 .- 1, linestyle = :dash, linewidth = lw)
    lines!(xs, ys .- i/6 .- 2, linestyle = :dot, linewidth = lw)
    lines!(xs, ys .- i/6 .- 3, linestyle = :dashdot, linewidth = lw)
    lines!(xs, ys .- i/6 .- 4, linestyle = :dashdotdot, linewidth = lw)
    lines!(xs, ys .- i/6 .- 5, linestyle = [0.5, 1.0, 1.5, 2.5], linewidth = lw)
end

f
```
\end{examplefigure}# series

{{doc series}}

## Examples

### Matrix

\begin{examplefigure}{}
```julia
using CairoMakie
CairoMakie.activate!() # hide
Makie.inline!(true) # hide

data = cumsum(randn(4, 101), dims = 2)

fig, ax, sp = series(data, labels=["label $i" for i in 1:4])
axislegend(ax)
fig
```
\end{examplefigure}

### Vector of vectors

\begin{examplefigure}{}
```julia
pointvectors = [Point2f.(1:100, cumsum(randn(100))) for i in 1:4]

series(pointvectors, markersize=5, color=:Set1)
```
\end{examplefigure}

### Vector and matrix

\begin{examplefigure}{}
```julia
data = cumsum(randn(4, 101), dims = 2)

series(0:0.1:10, data, solid_color=:black)
```
\end{examplefigure}
# stem

{{doc stem}}

### Examples

\begin{examplefigure}{}
```julia
using CairoMakie
CairoMakie.activate!() # hide
Makie.inline!(true) # hide

f = Figure()
Axis(f[1, 1])

xs = LinRange(0, 4pi, 30)

stem!(xs, sin.(xs))

f
```
\end{examplefigure}

\begin{examplefigure}{}
```julia
using CairoMakie
CairoMakie.activate!() # hide
Makie.inline!(true) # hide

f = Figure()
Axis(f[1, 1])

xs = LinRange(0, 4pi, 30)

stem!(xs, sin,
    offset = 0.5, trunkcolor = :blue, marker = :rect,
    stemcolor = :red, color = :orange,
    markersize = 15, strokecolor = :red, strokewidth = 3,
    trunklinestyle = :dash, stemlinestyle = :dashdot)

f
```
\end{examplefigure}

\begin{examplefigure}{}
```julia
using CairoMakie
CairoMakie.activate!() # hide
Makie.inline!(true) # hide

f = Figure()
Axis(f[1, 1])

xs = LinRange(0, 4pi, 30)

stem!(xs, sin.(xs),
    offset = LinRange(-0.5, 0.5, 30),
    color = LinRange(0, 1, 30), colorrange = (0, 0.5),
    trunkcolor = LinRange(0, 1, 30), trunkwidth = 5)

f
```
\end{examplefigure}

\begin{examplefigure}{}
```julia
using GLMakie
GLMakie.activate!() # hide
Makie.inline!(true) # hide

f = Figure()

xs = LinRange(0, 4pi, 30)

stem(f[1, 1], 0.5xs, 2 .* sin.(xs), 2 .* cos.(xs),
    offset = Point3f.(0.5xs, sin.(xs), cos.(xs)),
    stemcolor = LinRange(0, 1, 30), stemcolormap = :Spectral, stemcolorrange = (0, 0.5))

f
```
\end{examplefigure}
# boxplot

{{doc boxplot}}

### Examples

\begin{examplefigure}{}
```julia
using CairoMakie
CairoMakie.activate!() # hide
Makie.inline!(true) # hide

xs = rand(1:3, 1000)
ys = randn(1000)

boxplot(xs, ys)
```
\end{examplefigure}

\begin{examplefigure}{}
```julia
using CairoMakie
CairoMakie.activate!() # hide
Makie.inline!(true) # hide

xs = rand(1:3, 1000)
ys = randn(1000)
dodge = rand(1:2, 1000)

boxplot(xs, ys, dodge = dodge, show_notch = true, color = dodge)
```
\end{examplefigure}

Colors are customizable. The `color` attribute refers to the color of the boxes, whereas
`outliercolor` refers to the color of the outliers. If not scalars (e.g. `:red`), these attributes
must have the length of the data. If `outliercolor` is not provided, outliers will have the
same color as their box, as shown above.

!!! note
    For all indices corresponding to points within the same box, `color` (but not `outliercolor`)
    must have the same value.

\begin{examplefigure}{}
```julia
using CairoMakie
CairoMakie.activate!() # hide
Makie.inline!(true) # hide

xs = rand(1:3, 1000)
ys = randn(1000)
dodge = rand(1:2, 1000)

boxplot(xs, ys, dodge = dodge, show_notch = true, color = map(d->d==1 ? :blue : :red, dodge) , outliercolor = rand([:red, :green, :blue, :black, :yellow], 1000))
```
\end{examplefigure}
# hist

{{doc hist}}

### Examples

\begin{examplefigure}{}
```julia
using GLMakie
GLMakie.activate!() # hide
Makie.inline!(true) # hide

data = randn(1000)

f = Figure()
hist(f[1, 1], data, bins = 10)
hist(f[1, 2], data, bins = 20, color = :red, strokewidth = 1, strokecolor = :black)
hist(f[2, 1], data, bins = [-5, -2, -1, 0, 1, 2, 5], color = :gray)
hist(f[2, 2], data, normalization = :pdf)
f
```
\end{examplefigure}

#### Histogram with labels

You can use all the same arguments as \apilink{barplot}:
\begin{examplefigure}{}
```julia
using CairoMakie
CairoMakie.activate!() # hide

data = randn(1000)

hist(data, normalization = :pdf, bar_labels = :values,
     label_formatter=x-> round(x, digits=2), label_size = 15,
     strokewidth = 0.5, strokecolor = (:black, 0.5), color = :values)
```
\end{examplefigure}

#### Moving histograms

With `scale_to`, and `offset`, one can put multiple histograms into the same plot.
Note, that offset automatically sets fillto, to move the whole barplot.
Also, one can use a negative `scale_to` amount to flip the histogram.

\begin{examplefigure}{}
```julia
using CairoMakie
CairoMakie.activate!() # hide

fig = Figure()
ax = Axis(fig[1, 1])
for i in 1:5
     hist!(ax, randn(1000), scale_to=-0.6, offset=i, direction=:x)
end
fig
```
\end{examplefigure}
# surface

{{doc surface}}

## Attributes

### Generic

- `visible::Bool = true` sets whether the plot will be rendered or not.
- `overdraw::Bool = false` sets whether the plot will draw over other plots. This specifically means ignoring depth checks in GL backends.
- `transparency::Bool = false` adjusts how the plot deals with transparency. In GLMakie `transparency = true` results in using Order Independent Transparency.
- `fxaa::Bool = true` adjusts whether the plot is rendered with fxaa (anti-aliasing).
- `inspectable::Bool = true` sets whether this plot should be seen by `DataInspector`.
- `depth_shift::Float32 = 0f0` adjusts the depth value of a plot after all other transformations, i.e. in clip space, where `0 <= depth <= 1`. This only applies to GLMakie and WGLMakie and can be used to adjust render order (like a tunable overdraw). 
- `model::Makie.Mat4f` sets a model matrix for the plot. This replaces adjustments made with `translate!`, `rotate!` and `scale!`.
- `colormap::Union{Symbol, Vector{<:Colorant}} = :viridis` sets the colormap that is sampled for numeric `color`s.
- `colorrange::Tuple{<:Real, <:Real}` sets the values representing the start and end points of `colormap`.
- `nan_color::Union{Symbol, <:Colorant} = RGBAf(0,0,0,0)` sets a replacement color for `color = NaN`.

### Generic 3D

- `shading = true` enables lighting.
- `diffuse::Vec3f = Vec3f(0.4)` sets how strongly the red, green and blue channel react to diffuse (scattered) light. 
- `specular::Vec3f = Vec3f(0.2)` sets how strongly the object reflects light in the red, green and blue channels.
- `shininess::Real = 32.0` sets how sharp the reflection is.
- `ssao::Bool = false` adjusts whether the plot is rendered with ssao (screen space ambient occlusion). Note that this only makes sense in 3D plots and is only applicable with `fxaa = true`.

### Other

- `lowclip::Union{Nothing, Symbol, <:Colorant} = nothing` sets a color for any value below the colorrange.
- `highclip::Union{Nothing, Symbol, <:Colorant} = nothing` sets a color for any value above the colorrange.
- `invert_normals::Bool = false` inverts the normals generated for the surface. This can be useful to illuminate the other side of the surface.


## Examples

\begin{examplefigure}{}
```julia
using GLMakie
GLMakie.activate!() # hide
Makie.inline!(true) # hide

xs = LinRange(0, 10, 100)
ys = LinRange(0, 15, 100)
zs = [cos(x) * sin(y) for x in xs, y in ys]

surface(xs, ys, zs, axis=(type=Axis3,))
```
\end{examplefigure}

\begin{examplefigure}{}
```julia
using SparseArrays
using LinearAlgebra
using GLMakie
GLMakie.activate!() # hide
Makie.inline!(true) # hide

# This example was provided by Moritz Schauer (@mschauer).

#=
Define the precision matrix (inverse covariance matrix)
for the Gaussian noise matrix.  It approximately coincides
with the Laplacian of the 2d grid or the graph representing
the neighborhood relation of pixels in the picture,
https://en.wikipedia.org/wiki/Laplacian_matrix
=#
function gridlaplacian(m, n)
    S = sparse(0.0I, n*m, n*m)
    linear = LinearIndices((1:m, 1:n))
    for i in 1:m
        for j in 1:n
            for (i2, j2) in ((i + 1, j), (i, j + 1))
                if i2 <= m && j2 <= n
                    S[linear[i, j], linear[i2, j2]] -= 1
                    S[linear[i2, j2], linear[i, j]] -= 1
                    S[linear[i, j], linear[i, j]] += 1
                    S[linear[i2, j2], linear[i2, j2]] += 1
                end
            end
        end
    end
    return S
end

# d is used to denote the size of the data
d = 150

 # Sample centered Gaussian noise with the right correlation by the method
 # based on the Cholesky decomposition of the precision matrix
data = 0.1randn(d,d) + reshape(
        cholesky(gridlaplacian(d,d) + 0.003I) \ randn(d*d),
        d, d
)

surface(data; shading=false, colormap = :deep)
surface(data; shading=false, colormap = :deep)
```
\end{examplefigure}
# violin

{{doc violin}}

### Examples

\begin{examplefigure}{}
```julia
using CairoMakie
CairoMakie.activate!() # hide
Makie.inline!(true) # hide

xs = rand(1:3, 1000)
ys = randn(1000)

violin(xs, ys)
```
\end{examplefigure}

\begin{examplefigure}{}
```julia
using CairoMakie
CairoMakie.activate!() # hide
Makie.inline!(true) # hide

xs = rand(1:3, 1000)
ys = map(xs) do x
    return x == 1 ? randn() : x == 2 ? 0.5 * randn() : 5 * rand()
end

violin(xs, ys, datalimits = extrema)
```
\end{examplefigure}

\begin{examplefigure}{}
```julia
using CairoMakie
CairoMakie.activate!() # hide
Makie.inline!(true) # hide

N = 1000
xs = rand(1:3, N)
dodge = rand(1:2, N)
side = rand([:left, :right], N)
color = @. ifelse(side == :left, :orange, :teal)
ys = map(side) do s
    return s == :left ? randn() : rand()
end

violin(xs, ys, dodge = dodge, side = side, color = color)
```
\end{examplefigure}

\begin{examplefigure}{}
```julia
using CairoMakie
CairoMakie.activate!() # hide
Makie.inline!(true) # hide

N = 1000
xs = rand(1:3, N)
side = rand([:left, :right], N)
color = map(xs, side) do x, s
    colors = s == :left ? [:red, :orange, :yellow] : [:blue, :teal, :cyan]
    return colors[x]
end
ys = map(side) do s
    return s == :left ? randn() : rand()
end

violin(xs, ys, side = side, color = color)
```
\end{examplefigure}
# barplot

{{doc barplot}}

### Examples

\begin{examplefigure}{}
```julia
using CairoMakie
CairoMakie.activate!() # hide
Makie.inline!(true) # hide

f = Figure()
Axis(f[1, 1])

xs = 1:0.2:10
ys = 0.5 .* sin.(xs)

barplot!(xs, ys, color = :red, strokecolor = :black, strokewidth = 1)
barplot!(xs, ys .- 1, fillto = -1, color = xs, strokecolor = :black, strokewidth = 1)

f
```
\end{examplefigure}

\begin{examplefigure}{}
```julia
using CairoMakie
CairoMakie.activate!() # hide
Makie.inline!(true) # hide

xs = 1:0.2:10
ys = 0.5 .* sin.(xs)

barplot(xs, ys, gap = 0, color = :gray85, strokecolor = :black, strokewidth = 1)
```
\end{examplefigure}

\begin{examplefigure}{}
```julia
using CairoMakie
CairoMakie.activate!() # hide
Makie.inline!(true) # hide

tbl = (x = [1, 1, 1, 2, 2, 2, 3, 3, 3],
       height = 0.1:0.1:0.9,
       grp = [1, 2, 3, 1, 2, 3, 1, 2, 3],
       grp1 = [1, 2, 2, 1, 1, 2, 1, 1, 2],
       grp2 = [1, 1, 2, 1, 2, 1, 1, 2, 1]
       )

barplot(tbl.x, tbl.height,
        stack = tbl.grp,
        color = tbl.grp,
        axis = (xticks = (1:3, ["left", "middle", "right"]),
                title = "Stacked bars"),
        )
```
\end{examplefigure}

\begin{examplefigure}{}
```julia
barplot(tbl.x, tbl.height,
        dodge = tbl.grp,
        color = tbl.grp,
        axis = (xticks = (1:3, ["left", "middle", "right"]),
                title = "Dodged bars"),
        )
```
\end{examplefigure}

\begin{examplefigure}{}
```julia
barplot(tbl.x, tbl.height,
        dodge = tbl.grp1,
        stack = tbl.grp2,
        color = tbl.grp,
        axis = (xticks = (1:3, ["left", "middle", "right"]),
                title = "Dodged and stacked bars"),
        )
```
\end{examplefigure}

\begin{examplefigure}{}
```julia
colors = Makie.wong_colors()

# Figure and Axis
fig = Figure()
ax = Axis(fig[1,1], xticks = (1:3, ["left", "middle", "right"]),
        title = "Dodged bars with legend")

# Plot
barplot!(ax, tbl.x, tbl.height,
        dodge = tbl.grp,
        color = colors[tbl.grp])

# Legend
labels = ["group 1", "group 2", "group 3"]
elements = [PolyElement(polycolor = colors[i]) for i in 1:length(labels)]
title = "Groups"

Legend(fig[1,2], elements, labels, title)

fig
```
\end{examplefigure}

\begin{examplefigure}{}
```julia
barplot(
    tbl.x, tbl.height,
    dodge = tbl.grp,
    color = tbl.grp,
    bar_labels = :y,
    axis = (xticks = (1:3, ["left", "middle", "right"]),
            title = "Dodged bars horizontal with labels"),
    colormap = [:red, :green, :blue],
    color_over_background=:red,
    color_over_bar=:white,
    flip_labels_at=0.85,
    direction=:x,
)
```
\end{examplefigure}

\begin{examplefigure}{}
```julia
barplot(
    tbl.x, tbl.height,
    dodge = tbl.grp,
    color = tbl.grp,
    bar_labels = :y,
    axis = (xticks = (1:3, ["left", "middle", "right"]),
            title = "Dodged bars horizontal with labels"),
    colormap = [:red, :green, :blue],
    color_over_background=:red,
    color_over_bar=:white,
    flip_labels_at=0.85,
    direction=:x,
)
```
\end{examplefigure}

\begin{examplefigure}{}
```julia
barplot([-1, -0.5, 0.5, 1],
    bar_labels = :y,
    axis = (title="Fonts + flip_labels_at",),
    label_size = 20,
    flip_labels_at=(-0.8, 0.8),
    label_color=[:white, :green, :black, :white],
    label_formatter = x-> "Flip at $(x)?",
    label_offset = 10
)
```
\end{examplefigure}

\begin{examplefigure}{}
```julia
using CairoMakie

#Gantt data
gantt = (
    machine = [1,2,1,2],
    job = [1,1,2,3],
    task = [1,2,3,3],
    start = [1, 3, 3.5, 5],
    stop = [3, 4, 5, 6]
)

#Figure and axis
fig = Figure()
ax = Axis(
    fig[2,1],
    yticks = (1:2, ["A","B"]),
    ylabel = "Machine",
    xlabel = "Time"
)
xlims!(ax, 0, maximum(gantt.stop))

#Colors
colors = cgrad(:tab10)

#Plot bars
barplot!(
    gantt.machine,
    gantt.stop,
    fillto = gantt.start,
    direction = :x,
    color = colors[gantt.job],
    gap = 0.5
)

#Add labels
bar_labels = ["task #$i" for i in gantt.task]
text!(
    ["task #$i" for i in gantt.task],
    position = Point2f.(
        (gantt.start .+ gantt.stop) ./ 2,
        gantt.machine
    ),
    color = :white,
    align = (:center, :center)
)

#Add Legend
labels = ["job #$i" for i in unique(gantt.job)]
elements = [PolyElement(polycolor = colors[i]) for i in unique(gantt.job)]
Legend(fig[1,1], elements, labels, "Jobs", orientation=:horizontal, tellwidth = false, tellheight = true)

fig
```
\end{examplefigure}
# arrows

{{doc arrows}}

### Attributes

- `arrowhead = automatic`: Defines the marker (2D) or mesh (3D) that is used as
  the arrow head. The default for is `'▲'` in 2D and a cone mesh in 3D. For the
  latter the mesh should start at `Point3f(0)` and point in positive z-direction.
- `arrowtail = automatic`: Defines the mesh used to draw the arrow tail in 3D.
  It should start at `Point3f(0)` and extend in negative z-direction. The default
  is a cylinder. This has no effect on the 2D plot.
- `quality = 32`: Defines the number of angle subdivisions used when generating
  the arrow head and tail meshes. Consider lowering this if you have performance
  issues. Only applies to 3D plots.
- `linecolor = :black`: Sets the color used for the arrow tail which is
  represented by a line in 2D.
- `arrowcolor = linecolor`: Sets the color of the arrow head.
- `arrowsize = automatic`: Scales the size of the arrow head. This defaults to
  `0.3` in the 2D case and `Vec3f(0.2, 0.2, 0.3)` in the 3D case. For the latter
  the first two components scale the radius (in x/y direction) and the last scales
  the length of the cone. If the arrowsize is set to 1, the cone will have a
  diameter and length of 1.
- `linewidth = automatic`: Scales the width/diameter of the arrow tail.
  Defaults to `1` for 2D and `0.05` for the 3D case.
- `lengthscale = 1f0`: Scales the length of the arrow tail.
- `linestyle = nothing`: Sets the linestyle used in 2D. Does not apply to 3D
  plots.
- `normalize = false`: By default the lengths of the directions given to `arrows`
  are used to scale the length of the arrow tails. If this attribute is set to
  true the directions are normalized, skipping this scaling.
- `align = :origin`: Sets how arrows are positioned. By default arrows start at
  the given positions and extend along the given directions. If this attribute is
  set to `:head`, `:lineend`, `:tailend`, `:headstart` or `:center` the given
  positions will be between the head and tail of each arrow instead.

### Examples

\begin{examplefigure}{}
```julia
using CairoMakie
CairoMakie.activate!() # hide
Makie.inline!(true) # hide

f = Figure(resolution = (800, 800))
Axis(f[1, 1], backgroundcolor = "black")

xs = LinRange(0, 2pi, 20)
ys = LinRange(0, 3pi, 20)
us = [sin(x) * cos(y) for x in xs, y in ys]
vs = [-cos(x) * sin(y) for x in xs, y in ys]
strength = vec(sqrt.(us .^ 2 .+ vs .^ 2))

arrows!(xs, ys, us, vs, arrowsize = 10, lengthscale = 0.3,
    arrowcolor = strength, linecolor = strength)

f
```
\end{examplefigure}

\begin{examplefigure}{}
```julia
using GLMakie
GLMakie.activate!() # hide
Makie.inline!(true) # hide

ps = [Point3f(x, y, z) for x in -5:2:5 for y in -5:2:5 for z in -5:2:5]
ns = map(p -> 0.1 * Vec3f(p[2], p[3], p[1]), ps)
arrows(
    ps, ns, fxaa=true, # turn on anti-aliasing
    linecolor = :gray, arrowcolor = :black,
    linewidth = 0.1, arrowsize = Vec3f(0.3, 0.3, 0.4),
    align = :center, axis=(type=Axis3,)
)
```
\end{examplefigure}

\begin{examplefigure}{}
```julia
using GLMakie
GLMakie.activate!() # hide
Makie.inline!(true) # hide
using LinearAlgebra

ps = [Point3f(x, y, z) for x in -5:2:5 for y in -5:2:5 for z in -5:2:5]
ns = map(p -> 0.1 * Vec3f(p[2], p[3], p[1]), ps)
lengths = norm.(ns)
arrows(
    ps, ns, fxaa=true, # turn on anti-aliasing
    color=lengths,
    linewidth = 0.1, arrowsize = Vec3f(0.3, 0.3, 0.4),
    align = :center, axis=(type=Axis3,)
)
```
\end{examplefigure}
# hspan! and vspan!

{{doc hspan!}}
{{doc vspan!}}

These functions are not plot types / recipes and only work with `Axis`.

\begin{examplefigure}{}
```julia
using CairoMakie
CairoMakie.activate!() # hide

f = Figure()
ax = Axis(f[1, 1])

lines!(ax, 0..20, sin)
vspan!(ax, [0, 2pi, 4pi], [pi, 3pi, 5pi], color = (:red, 0.2))
hspan!(ax, -1.1, -0.9, color = (:blue, 0.2))

f
```
\end{examplefigure}
# scatterlines

{{doc scatterlines}}

## Examples

\begin{examplefigure}{svg = true}
```julia
using CairoMakie
CairoMakie.activate!() # hide
Makie.inline!(true) # hide

f = Figure()
Axis(f[1, 1])

xs = LinRange(0, 10, 20)
ys = 0.5 .* sin.(xs)

scatterlines!(xs, ys, color = :red)
scatterlines!(xs, ys .- 1, color = xs, markercolor = :red)
scatterlines!(xs, ys .- 2, markersize = LinRange(5, 30, 20))
scatterlines!(xs, ys .- 3, marker = :cross, strokewidth = 1,
    markersize = 20, color = :orange, strokecolor = :black)

f
```
\end{examplefigure}
# stairs

{{doc stairs}}

### Examples

\begin{examplefigure}{}
```julia
using CairoMakie
CairoMakie.activate!() # hide
Makie.inline!(true) # hide

f = Figure()

xs = LinRange(0, 4pi, 21)
ys = sin.(xs)

stairs(f[1, 1], xs, ys)
stairs(f[2, 1], xs, ys; step=:post, color=:blue, linestyle=:dash)
stairs(f[3, 1], xs, ys; step=:center, color=:red, linestyle=:dot)

f
```
\end{examplefigure}
# wireframe

{{doc wireframe}}

### Examples

\begin{examplefigure}{}
```julia
using GLMakie
GLMakie.activate!() # hide
Makie.inline!(true) # hide

x, y = collect(-8:0.5:8), collect(-8:0.5:8)
z = [sinc(√(X^2 + Y^2) / π) for X ∈ x, Y ∈ y]

wireframe(x, y, z, axis=(type=Axis3,), color=:black)
```
\end{examplefigure}
# heatmap

{{doc heatmap}}

## Attributes

### Generic

- `visible::Bool = true` sets whether the plot will be rendered or not.
- `overdraw::Bool = false` sets whether the plot will draw over other plots. This specifically means ignoring depth checks in GL backends.
- `transparency::Bool = false` adjusts how the plot deals with transparency. In GLMakie `transparency = true` results in using Order Independent Transparency.
- `fxaa::Bool = true` adjusts whether the plot is rendered with fxaa (anti-aliasing).
- `inspectable::Bool = true` sets whether this plot should be seen by `DataInspector`.
- `depth_shift::Float32 = 0f0` adjusts the depth value of a plot after all other transformations, i.e. in clip space, where `0 <= depth <= 1`. This only applies to GLMakie and WGLMakie and can be used to adjust render order (like a tunable overdraw).
- `model::Makie.Mat4f` sets a model matrix for the plot. This replaces adjustments made with `translate!`, `rotate!` and `scale!`.
- `color` is set by the plot.
- `colormap::Union{Symbol, Vector{<:Colorant}} = :viridis` sets the colormap that is sampled for numeric `color`s.
- `colorrange::Tuple{<:Real, <:Real}` sets the values representing the start and end points of `colormap`.
- `nan_color::Union{Symbol, <:Colorant} = RGBAf(0,0,0,0)` sets a replacement color for `color = NaN`.

### Other

- `lowclip::Union{Nothing, Symbol, <:Colorant} = nothing` sets a color for any value below the colorrange.
- `highclip::Union{Nothing, Symbol, <:Colorant} = nothing` sets a color for any value above the colorrange.
- `interpolate::Bool = false` sets whether colors should be interpolated.


## Examples

### Two vectors and a matrix

\begin{examplefigure}{}
```julia
using CairoMakie
CairoMakie.activate!() # hide
Makie.inline!(true) # hide


xs = range(0, 10, length = 25)
ys = range(0, 15, length = 25)
zs = [cos(x) * sin(y) for x in xs, y in ys]

heatmap(xs, ys, zs)
```
\end{examplefigure}

### Two ranges and a function

\begin{examplefigure}{name = "mandelbrot_heatmap"}
```julia
using CairoMakie
CairoMakie.activate!() # hide
Makie.inline!(true) # hide

function mandelbrot(x, y)
    z = c = x + y*im
    for i in 1:30.0; abs(z) > 2 && return i; z = z^2 + c; end; 0
end

heatmap(-2:0.1:1, -1.1:0.1:1.1, mandelbrot,
    colormap = Reverse(:deep))
```
\end{examplefigure}

### Three vectors

There must be no duplicate combinations of x and y, but it is allowed to leave out values.

\begin{examplefigure}{}
```julia
using CairoMakie
CairoMakie.activate!() # hide
Makie.inline!(true) # hide


xs = [1, 2, 3, 1, 2, 3, 1, 2, 3]
ys = [1, 1, 1, 2, 2, 2, 3, 3, 3]
zs = [1, 2, 3, 4, 5, 6, 7, 8, NaN]

heatmap(xs, ys, zs)
```
\end{examplefigure}

### Colors

Using the previous example, we can see how to change the color of our plot. You can find additional colors [here](https://makie.juliaplots.org/stable/documentation/colors/index.html#colors).
# Basic Tutorial

## Preface

Here is a quick tutorial to get you started with Makie!

We assume you have [Julia](https://julialang.org/) and `CairoMakie.jl` (or one of the other backends, i.e. `GLMakie.jl` or `WGLMakie.jl`) installed already.

This tutorial uses CairoMakie, but the code can be executed with any backend.
CairoMakie can output beautiful static vector graphics, but it doesn't have the native ability to open interactive windows.

To see the output of plotting commands when using CairoMakie, we recommend you either use an IDE which supports png or svg output, such as VSCode, Atom/Juno, Jupyter, Pluto, etc., or try using a viewer package such as [ElectronDisplay.jl](https://github.com/queryverse/ElectronDisplay.jl), or alternatively save your plots to files directly.
The Julia REPL by itself does not have the ability to show the plots.

GLMakie can open interactive plot windows, also from the Julia REPL, or alternatively display bitmaps inline if `Makie.inline!(true)` is called
and if it is supported by the environment.

WGLMakie shows interactive plots in environments that support interactive html displays, such as VSCode, Atom/Juno, Jupyter, Pluto, etc.

For more information, have a look at \myreflink{Backends & Output}.

Ok, now that this is out of the way, let's get started!

## First plot

First, we import CairoMakie.

```julia:setup
using CairoMakie
CairoMakie.activate!() # hide
Makie.inline!(true) # hide
nothing # hide
```

Makie has many different plotting functions, one of the most common ones is \myreflink{lines}.
You can just call such a function and your plot will appear if your coding environment can show png or svg files.

!!! note
    Objects such as [Figure](\reflink{Figures}), `FigureAxisPlot` or `Scene` are usually displayed whenever they are returned in global scope (e.g. in the REPL).
    To display such objects from within a local scope, like from within a function, you can directly call `display(figure)`, for example.

\begin{examplefigure}{svg = true}
```julia

x = range(0, 10, length=100)
y = sin.(x)
lines(x, y)
```
\end{examplefigure}

Another common function is \myreflink{scatter}.

\begin{examplefigure}{svg = true}
```julia

using CairoMakie

x = range(0, 10, length=100)
y = sin.(x)
scatter(x, y)
```
\end{examplefigure}

## Multiple plots

Every plotting function has a version with and one without `!`.
For example, there's `scatter` and `scatter!`, `lines` and `lines!`, etc.
The functions without a `!` always create a new axis with a plot inside, while the functions with `!` plot into an already existing axis.

Here's how you could plot two lines on top of each other.

\begin{examplefigure}{svg = true}
```julia

using CairoMakie

x = range(0, 10, length=100)
y1 = sin.(x)
y2 = cos.(x)

lines(x, y1)
lines!(x, y2)
current_figure()
```
\end{examplefigure}

The second `lines!` call plots into the axis created by the first `lines` call.
If you don't specify an axis to plot into, it's as if you had called `lines!(current_axis(), ...)`.

The call to `current_figure` is necessary here, because functions with `!` return only the newly created plot object, but this alone does not cause the figure to display when returned.

## Attributes

Every plotting function has attributes which you can set through keyword arguments.
The lines in the previous example both have the same default color, which we can change easily.

\begin{examplefigure}{svg = true}
```julia

using CairoMakie

x = range(0, 10, length=100)
y1 = sin.(x)
y2 = cos.(x)

lines(x, y1, color = :red)
lines!(x, y2, color = :blue)
current_figure()
```
\end{examplefigure}

Other plotting functions have different attributes.
The function `scatter`, for example, does not only have the `color` attribute, but also a `markersize` attribute.

\begin{examplefigure}{svg = true}
```julia

using CairoMakie

x = range(0, 10, length=100)
y1 = sin.(x)
y2 = cos.(x)

scatter(x, y1, color = :red, markersize = 5)
scatter!(x, y2, color = :blue, markersize = 10)
current_figure()
```
\end{examplefigure}

If you save the plot object returned from a call like `scatter!`, you can also manipulate its attributes later with the syntax `plot.attribute = new_value`.

\begin{examplefigure}{svg = true}
```julia

using CairoMakie

x = range(0, 10, length=100)
y1 = sin.(x)
y2 = cos.(x)

scatter(x, y1, color = :red, markersize = 5)
sc = scatter!(x, y2, color = :blue, markersize = 10)
sc.color = :green
sc.markersize = 20
current_figure()
```
\end{examplefigure}

## Array attributes

A lot of attributes can be set to either a single value or an array with as many elements as there are data points.
For example, it is usually much more performant to draw many points with one scatter object, than to create many scatter objects with one point each.

Here are the two scatter plots again, but one has varying markersize, and the other varying color.

\begin{examplefigure}{svg = true}
```julia

using CairoMakie

x = range(0, 10, length=100)
y1 = sin.(x)
y2 = cos.(x)

scatter(x, y1, color = :red, markersize = range(5, 15, length=100))
sc = scatter!(x, y2, color = range(0, 1, length=100), colormap = :thermal)

current_figure()
```
\end{examplefigure}

Note that the color array does not actually contain colors, rather the numerical values are mapped to the plot's `colormap`.
There are many different colormaps to choose from, take a look on the \myreflink{Colors} page.

The values are mapped to colors via the `colorrange` attribute, which by default goes from the minimum to the maximum color value, but we can also limit or expand the range manually.
For example, we can constrain the previous scatter plot's color range to (0.25, 0.75), which will clip the colors at the bottom and the top quarters.

\begin{examplefigure}{svg = true}
```julia

sc.colorrange = (0.25, 0.75)

current_figure()
```
\end{examplefigure}

Of course you can also use an array of colors directly, in which case the `colorrange` is ignored:

\begin{examplefigure}{svg = true}
```julia

using CairoMakie

x = range(0, 10, length=100)
y = sin.(x)

colors = repeat([:crimson, :dodgerblue, :slateblue1, :sienna1, :orchid1], 20)

scatter(x, y, color = colors, markersize = 20)
```
\end{examplefigure}

## Simple legend

If you add label attributes to your plots, you can call the `axislegend` function to add a legend with all labeled plots to the current axis.

\begin{examplefigure}{svg = true}
```julia

using CairoMakie

x = range(0, 10, length=100)
y1 = sin.(x)
y2 = cos.(x)

lines(x, y1, color = :red, label = "sin")
lines!(x, y2, color = :blue, label = "cos")
axislegend()
current_figure()
```
\end{examplefigure}

## Subplots

Makie uses a powerful layout system under the hood, which allows you to create very complex figures with many subplots.
For the easiest way to do this, we need a [Figure](\reflink{Figures}) object.
So far, we haven't seen this explicitly, it was created in the background in the first plotting function call.

We can also create a [Figure](\reflink{Figures}) directly and then continue working with it.
We can make subplots by giving the location of the subplot in our layout grid as the first argument to our plotting function.
The basic syntax for specifying the location in a figure is `fig[row, col]`.

\begin{examplefigure}{svg = true}
```julia

using CairoMakie

x = LinRange(0, 10, 100)
y = sin.(x)

fig = Figure()
lines(fig[1, 1], x, y, color = :red)
lines(fig[1, 2], x, y, color = :blue)
lines(fig[2, 1:2], x, y, color = :green)

fig
```
\end{examplefigure}

Each `lines` call creates a new axis in the position given as the first argument, that's why we use `lines` and not `lines!` here.

## Constructing axes manually

Like [Figure](\reflink{Figures})s, we can also create axes manually.
This is useful if we want to prepare an empty axis to then plot into it later.

The default 2D axis that we have created implicitly so far is called \myreflink{Axis} and can also be created in a specific position in the figure by passing that position as the first argument.

For example, we can create a figure with three axes.

\begin{examplefigure}{svg = true}
```julia

using CairoMakie

fig = Figure()
ax1 = Axis(fig[1, 1])
ax2 = Axis(fig[1, 2])
ax3 = Axis(fig[2, 1:2])
fig
```
\end{examplefigure}

And then we can continue to plot into these empty axes.

\begin{examplefigure}{svg = true}
```julia

lines!(ax1, 0..10, sin)
lines!(ax2, 0..10, cos)
lines!(ax3, 0..10, sqrt)
fig
```
\end{examplefigure}

Note, the notation `0..10` above creates a closed interval from `0` to `10` (see [`IntervalSets.jl`](https://github.com/JuliaMath/IntervalSets.jl) for further details).

Axes also have many attributes that you can set, for example to give them a title, or labels.

\begin{examplefigure}{svg = true}
```julia

ax1.title = "sin"
ax2.title = "cos"
ax3.title = "sqrt"

ax1.ylabel = "amplitude"
ax3.ylabel = "amplitude"
ax3.xlabel = "time"
fig
```
\end{examplefigure}

## Legend and Colorbar

We have seen two `Layoutables` so far, the \myreflink{Axis} and the \myreflink{Legend} which was created by the function `axislegend`.
All `Layoutable`s can be placed into the layout of a figure at arbitrary positions, which makes it easy to assemble complex figures.

In the same way as with the \myreflink{Axis} before, you can also create a \myreflink{Legend} manually and then place it freely, wherever you want, in the figure.
There are multiple ways to create \myreflink{Legend}s, for one of them you pass one vector of plot objects and one vector of label strings.

You can see here that we can deconstruct the return value from the two `lines` calls into one newly created axis and one plot object each.
We can then feed the plot objects to the legend constructor.
We place the legend in the second column and across both rows, which centers it nicely next to the two axes.

\begin{examplefigure}{svg = true}
```julia

using CairoMakie

fig = Figure()
ax1, l1 = lines(fig[1, 1], 0..10, sin, color = :red)
ax2, l2 = lines(fig[2, 1], 0..10, cos, color = :blue)
Legend(fig[1:2, 2], [l1, l2], ["sin", "cos"])
fig
```
\end{examplefigure}

The \myreflink{Colorbar} works in a very similar way.
We just need to pass a position in the figure to it, and one plot object.
In this example, we use a `heatmap`.

You can see here that we split the return value of `heatmap` into three parts: the newly created figure, the axis and the heatmap plot object.
This is useful as we can then continue with the figure `fig` and the heatmap `hm` which we need for the colorbar.

\begin{examplefigure}{svg = true}
```julia

using CairoMakie

fig, ax, hm = heatmap(randn(20, 20))
Colorbar(fig[1, 2], hm)
fig
```
\end{examplefigure}

The previous short syntax is basically equivalent to this longer, manual version.
You can switch between those workflows however you please.

\begin{examplefigure}{svg = true}
```julia

using CairoMakie

fig = Figure()
ax = Axis(fig[1, 1])
hm = heatmap!(ax, randn(20, 20))
Colorbar(fig[1, 2], hm)
fig
```
\end{examplefigure}

## Passing attributes to Figure and Axis

For one-off plots, it can be convenient to set axis or figure settings directly with the plotting command.
You can do this using the plotting functions without the `!` suffix, like `lines` or `scatter`, because these always create a new axis and also create a new figure if they are not plotting onto an existing one. This is explained further under \myreflink{Plot Method Signatures}.

You can pass axis attributes under the keyword `axis` and figure attributes under the keyword `figure`.

\begin{examplefigure}{svg = true}
```julia

using CairoMakie

heatmap(randn(20, 20),
    figure = (backgroundcolor = :pink,),
    axis = (aspect = 1, xlabel = "x axis", ylabel = "y axis")
)
```
\end{examplefigure}

If you set only one attribute, be careful to do `axis = (key = value,)` (note the trailing comma), otherwise you're not creating a `NamedTuple` but a local variable `key`.

## Next steps

We've only looked at a small subset of Makie's functionality here.

You can read about the different available plotting functions with examples in the Plotting Functions section.

If you want to learn about making complex figures with nested sublayouts, have a look at the \myreflink{Layout Tutorial} section.

If you're interested in creating interactive visualizations that use Makie's special `Observables` workflow, this is explained in more detail in the \myreflink{Observables & Interaction} section.

If you want to create animated movies, you can find more information in the \myreflink{Animations} section.
@def reeval = false

# Layout Tutorial

In this tutorial, you will learn how to create a complex figure using Makie's layout tools.

Let's say that we want to create the following figure:

\outputimage{final_result.png}

How do we approach this task?

In the following sections, we'll go over the process step by step.
We're not always going to use the shortest possible syntax, as the main goal is to get a better understanding of the logic and the available options.

## Basic layout plan

When building figures, you always think in terms of rectangular boxes. We want to find the biggest boxes that enclose meaningful groups of content, and then we realize those boxes either using `GridLayout` or by placing content objects there.

If we look at our target figure, we can imagine one box around each of the labelled areas A, B, C and D. But A and C are not in one row, neither are B and D. This means that we don't use a 2x2 GridLayout, but have to be a little more creative.

We could say that A and B are in one column, and C and D are in one column. We can have different row heights for both groups by making one big nested `GridLayout` within the second column, in which we place C and D. This way the rows of column 2 are decoupled from column 1.

Ok, let's create the figure first with a gray backgroundcolor, and a predefined font:

\begin{examplefigure}{}
```julia
using CairoMakie
using FileIO
CairoMakie.activate!() # hide

noto_sans = assetpath("fonts", "NotoSans-Regular.ttf")
noto_sans_bold = assetpath("fonts", "NotoSans-Bold.ttf")

f = Figure(backgroundcolor = RGBf(0.98, 0.98, 0.98),
    resolution = (1000, 700), font = noto_sans)
```
\end{examplefigure}

## Setting up GridLayouts

Now, let's make the four nested GridLayouts that are going to hold the objects of A, B, C and D. There's also the layout that holds C and D together, so the rows are separate from A and B. We are not going to see anything yet as we have no visible content, but that will come soon.

!!! note
    It's not strictly necessary to first create separate `GridLayout`s, then use them to place objects in the figure. You can also implicitly create nested grids using multiple indexing, for example like `Axis(f[1, 2:3][4:5, 6])`. This is further explained in \myreflink{GridPositions and GridSubpositions}. But if you want to manipulate your nested grids afterwards, for example to change column sizes or row gaps, it's easier if you have them stored in variables already.

```julia:grids
ga = f[1, 1] = GridLayout()
gb = f[2, 1] = GridLayout()
gcd = f[1:2, 2] = GridLayout()
gc = gcd[1, 1] = GridLayout()
gd = gcd[2, 1] = GridLayout()
```

## Panel A

Now we can start placing objects into the figure. We start with A.

There are three axes and a legend. We can place the axes first, link them appropriately, and plot the first data into them.

\begin{examplefigure}{}
```julia
axtop = Axis(ga[1, 1])
axmain = Axis(ga[2, 1], xlabel = "before", ylabel = "after")
axright = Axis(ga[2, 2])

linkyaxes!(axmain, axright)
linkxaxes!(axmain, axtop)

labels = ["treatment", "placebo", "control"]
data = randn(3, 100, 2) .+ [1, 3, 5]

for (label, col) in zip(labels, eachslice(data, dims = 1))
    scatter!(axmain, col, label = label)
    density!(axtop, col[:, 1])
    density!(axright, col[:, 2], direction = :y)
end

f
```
\end{examplefigure}

There's a small gap between the density plots and their axes, which we can remove by fixing one side of the limits.

\begin{examplefigure}{}
```julia
ylims!(axtop, low = 0)
xlims!(axright, low = 0)

f
```
\end{examplefigure}

### Legend

We have set the `label` attribute in the scatter call so it's easier to construct the legend. We can just pass `axmain` as the second argument to `Legend`.

\begin{examplefigure}{}
```julia
leg = Legend(ga[1, 2], axmain)

f
```
\end{examplefigure}

### Legend Tweaks

There are a couple things we want to change. There are unnecessary decorations for the side axes, which we are going to hide.

Also, the top axis does not have the same height as the legend. That's because a legend is usually used on the right of an `Axis` and is therefore preset with `tellheight = false`. We set this attribute to `true` so the row in which the legend sits can contract to its known size.

\begin{examplefigure}{}
```julia
hidedecorations!(axtop, grid = false)
hidedecorations!(axright, grid = false)
leg.tellheight = true

f
```
\end{examplefigure}

The axes are still a bit too far apart, so we reduce column and row gaps.

\begin{examplefigure}{}
```julia
colgap!(ga, 10)
rowgap!(ga, 10)

f
```
\end{examplefigure}

We can make a title by placing a label across the top two elements.

\begin{examplefigure}{}
```julia
Label(ga[1, 1:2, Top()], "Stimulus ratings", valign = :bottom,
    padding = (0, 0, 5, 0))

f
```
\end{examplefigure}

## Panel B

Let's move to B. We have two axes stacked on top of each other, and a colorbar alongside them. This time, we create the axes by just plotting into the right `GridLayout` slots. This can be more convenient than creating an `Axis` first.

\begin{examplefigure}{}
```julia
xs = LinRange(0.5, 6, 50)
ys = LinRange(0.5, 6, 50)
data1 = [sin(x^1.5) * cos(y^0.5) for x in xs, y in ys] .+ 0.1 .* randn.()
data2 = [sin(x^0.8) * cos(y^1.5) for x in xs, y in ys] .+ 0.1 .* randn.()

ax1, hm = contourf(gb[1, 1], xs, ys, data1,
    levels = 6)
ax1.title = "Histological analysis"
contour!(xs, ys, data1, levels = 5, color = :black)
hidexdecorations!(ax1)

_, hm2 = contourf(gb[2, 1], xs, ys, data2,
    levels = 6)
contour!(xs, ys, data2, levels = 5, color = :black)

f
```
\end{examplefigure}

### Colorbar

Now we need a colorbar.
Because we haven't set specific edges for the two contour plots, just how many levels there are, we can make a colorbar using one of the contour plots and then label the bins in there from one to six.

\begin{examplefigure}{}
```julia
cb = Colorbar(gb[1:2, 2], hm, label = "cell group")
low, high = extrema(data1)
edges = range(low, high, length = 7)
centers = (edges[1:6] .+ edges[2:7]) .* 0.5
cb.ticks = (centers, string.(1:6))

f
```
\end{examplefigure}

#### Mixed alignmode

The right edge of the colorbar is currently aligned with the right edge of the upper density plot.
This can later cause a bit of a gap between the density plot and content on the right.

In order to improve this, we can pull the colorbar labels into its layout cell using the `Mixed` alignmode. The keyword `right = 0` means that the right side of the colorbar should pull its protrusion content inward with an additional padding of `0`.

\begin{examplefigure}{}
```julia
cb.alignmode = Mixed(right = 0)

f
```
\end{examplefigure}

As in A, the axes are a bit too far apart.

\begin{examplefigure}{}
```julia
colgap!(gb, 10)
rowgap!(gb, 10)

f
```
\end{examplefigure}

## Panel C

Now, we move on to panel C. This is just an `Axis3` with a colorbar on the side.

\begin{examplefigure}{}
```julia
brain = load(assetpath("brain.stl"))

Axis3(gc[1, 1], title = "Brain activation")
m = mesh!(
    brain,
    color = [tri[1][2] for tri in brain for i in 1:3],
    colormap = Reverse(:magma),
)
Colorbar(gc[1, 2], m, label = "BOLD level")

f
```
\end{examplefigure}

Note that the z label overlaps the plot to the left a little bit. `Axis3` can't have automatic protrusions because the label positions change with the projection and the cell size of the axis, which is different from the 2D `Axis`.

You can set the attribute `ax3.protrusions` to a tuple of four values (left, right, bottom, top) but in this case we just continue plotting until we have all objects that we want, before we look if small tweaks like that are necessary.

## Panel D

We move on to Panel D, which has a grid of 3x2 axes.

\begin{examplefigure}{}
```julia
axs = [Axis(gd[row, col]) for row in 1:3, col in 1:2]
hidedecorations!.(axs, grid = false, label = false)

for row in 1:3, col in 1:2
    xrange = col == 1 ? (0:0.1:6pi) : (0:0.1:10pi)

    eeg = [sum(sin(pi * rand() + k * x) / k for k in 1:10)
        for x in xrange] .+ 0.1 .* randn.()

    lines!(axs[row, col], eeg, color = (:black, 0.5))
end

axs[3, 1].xlabel = "Day 1"
axs[3, 2].xlabel = "Day 2"

f
```
\end{examplefigure}

We can make a little title for the six axes by placing a `Label` in the top protrusion of row 1 and across both columns.

\begin{examplefigure}{}
```julia
Label(gd[1, :, Top()], "EEG traces", valign = :bottom,
    padding = (0, 0, 5, 0))

f
```
\end{examplefigure}

Again, we bring the subplots closer together by reducing gap sizes.

\begin{examplefigure}{}
```julia
rowgap!(gd, 10)
colgap!(gd, 10)

f
```
\end{examplefigure}

### EEG labels

Now, we add three boxes on the side with labels in them. In this case, we just place them in another column to the right.

\begin{examplefigure}{}
```julia
for (i, label) in enumerate(["sleep", "awake", "test"])
    Box(gd[i, 3], color = :gray90)
    Label(gd[i, 3], label, rotation = pi/2, tellheight = false)
end

f
```
\end{examplefigure}

The boxes are in the correct positions, but we still need to remove the column gap.

\begin{examplefigure}{}
```julia
colgap!(gd, 2, 0)

f
```
\end{examplefigure}

### Scaling axes relatively

The fake eeg data we have created has more datapoints on day 1 than day 2.
We want to scale the axes so that they both have the same zoom level.
We can do this by setting the column widths to `Auto(x)` where x is a number proportional to the number of data points of the axis.
This way, both will have the same relative scaling.

\begin{examplefigure}{}
```julia
n_day_1 = length(0:0.1:6pi)
n_day_2 = length(0:0.1:10pi)

colsize!(gd, 1, Auto(n_day_1))
colsize!(gd, 2, Auto(n_day_2))

f
```
\end{examplefigure}

## Subplot labels

Now, we can add the subplot labels. We already have our four `GridLayout` objects that enclose each panel's content, so the easiest way is to create `Label`s in the top left protrusion of these layouts.
That will leave all other alignments intact, because we're not creating any new columns or rows. The labels belong to the gaps between the layouts instead.

\begin{examplefigure}{}
```julia
for (label, layout) in zip(["A", "B", "C", "D"], [ga, gb, gc, gd])
    Label(layout[1, 1, TopLeft()], label,
        textsize = 26,
        font = noto_sans_bold,
        padding = (0, 5, 5, 0),
        halign = :right)
end

f
```
\end{examplefigure}

## Final tweaks

This looks pretty good already, but the first column of the layout is a bit too wide.
We can reduce the column width by setting it to `Auto` with a number smaller than 1, for example.
This gives the column a smaller weight when distributing widths between all columns with `Auto` sizes.

You can also use `Relative` or `Fixed` but they are not as flexible if you add more things later, so I prefer using `Auto`.

\begin{examplefigure}{}
```julia
colsize!(f.layout, 1, Auto(0.5))

f
```
\end{examplefigure}

The EEG traces are currently as high as the brain axis, let's increase the size of the row with the panel C layout a bit so it has more space.

And that is the final result:

\begin{examplefigure}{name = "final_result"}
```julia
rowsize!(gcd, 1, Auto(1.5))

f
```
\end{examplefigure}
# Aspect ratio and size control tutorial

A very common problem in plotting is dealing with aspect ratios and other ways to precisely control figures.

For example, many plots need square axes.
If you have looked at the documentation of `Axis`, you might know that it has an `aspect` attribute that can control the aspect ratio of the axis box.
This aspect is not concerned with what the data limits are, it's just about the relative visual length of the axes.

Let's look at one common example, a square axis with a colorbar next to it:


\begin{examplefigure}{svg = true}
```julia
using CairoMakie
CairoMakie.activate!() # hide

set_theme!(backgroundcolor = :gray90)

f = Figure(resolution = (800, 500))
ax = Axis(f[1, 1], aspect = 1)
Colorbar(f[1, 2])
f
```
\end{examplefigure}


As you can see, the axis is square, but there's also a large gap between it and the colorbar.
Why is that?

We can visualize the reason by adding a Box to the same cell where the axis is:


\begin{examplefigure}{svg = true, name = "aspect_tutorial_example"}
```julia
Box(f[1, 1], color = (:red, 0.2), strokewidth = 0)
f
```
\end{examplefigure}


The red area of the box extends out into the whitespace left by the Axis.
This demonstrates what the `aspect` keyword is actually doing.
It reduces the size of the Axis, such that the chosen aspect ratio is achieved.
It doesn't tell the layout that the Axis lives in "please make this cell adhere to this aspect ratio".
As far as the layout is concerned, the Axis has an undefined size and its layout cell can therefore have any size that the layout deems correct, based on all other content of the layout and the figure size.

Therefore, using `aspect` will always cause gaps, unless the layout cell where the Axis lives happens to have exactly the correct aspect ratio by chance.
This means `aspect` should only be used if the whitespace caused by it does not matter too much.

For all other cases, there is a different approach.

We want to force the layout to keep the axis cell at a specific aspect ratio.
Therefore, we have to manipulate the layout itself, not the axis.

By default, each GridLayout row and column has a size of `Auto()`.
This means that the size can depend on fixed-size content if there is any, otherwise it expands to fill the available space.
If we want to force a cell to have an aspect ratio, we need to set either its respective row or column size to `Aspect`.

Let's try the example from above again, but this time we force the column of the Axis to have an aspect ratio of 1.0 relative to the row of the Axis, which is row 1.


\begin{examplefigure}{svg = true}
```julia
f = Figure(resolution = (800, 500))
ax = Axis(f[1, 1])
Colorbar(f[1, 2])
colsize!(f.layout, 1, Aspect(1, 1.0))
f
```
\end{examplefigure}


As you can see, this time the colorbar sticks close to the axis, there is no unnecessary whitespace between them.
We can visualize the effect of `Aspect` again with a red box, that shows us the extent of the layout cell:


\begin{examplefigure}{svg = true}
```julia
# hide
Box(f[1, 1], color = (:red, 0.2), strokewidth = 0)
f
```
\end{examplefigure}


So this time the layout cell itself is square, therefore the Axis that fills it is also square.
Let me just demonstrate that we can play the same game again and give the Axis an `aspect` that is different from the square one that the layout cell has.
This will again cause unnecessary whitespace:


\begin{examplefigure}{svg = true}
```julia
ax.aspect = 0.5
f
```
\end{examplefigure}


And now we change the column aspect again, to remove this gap:


\begin{examplefigure}{svg = true}
```julia
colsize!(f.layout, 1, Aspect(1, 0.5))
f
```
\end{examplefigure}


Let's return to our previous state with a square axis:


\begin{examplefigure}{svg = true}
```julia
# hide
f = Figure(resolution = (800, 500))
ax = Axis(f[1, 1])
Colorbar(f[1, 2])
colsize!(f.layout, 1, Aspect(1, 1.0))
f
```
\end{examplefigure}


Now you might think that there is no whitespace anymore between Axis and Colorbar, but there is a lot of it to the left and the right.
Why can the layout not fix this problem for us?

Well, in Makie, the layout has to operate within the confines of the figure size that we have set.
It cannot just decrease the figure size if there's too little content.
This is because lots of times, figures are created to fit the sizing rules of some academic journal exactly, therefore the content you plot is not allowed to mess with the figure size.

So what we have done in our example is introducing constraints to the sizes of objects in our layout, such that it's impossible to fill all the space that is theoretically available.
If you think about it, it's impossible to fill this Figure with a square axis and a thin colorbar while filling the rectangular space.
We need a smaller figure!

But how small should it be exactly?
It would be quite difficult to eyeball this, but thankfully there's a function for this exact purpose.
By calling `resize_to_layout!`, we can adjust the figure size to the size that the layout needs for all its content.

Let's try it out:


\begin{examplefigure}{svg = true}
```julia
resize_to_layout!(f)
f
```
\end{examplefigure}


As you can see, the whitespace at the sides has been trimmed.
(If the scaling looks smaller or bigger, that is just because of the display on this site, not the underlying figure size).

This technique is useful for all kinds of situations where the content should decide the figure size, and not the other way around.

For example, let's say we have a facet plot with 25 square axes which are all of size 150 by 150.
We can just make these axes with fixed widths and heights.
The `Auto` sized columns and rows of the default layout pick up these measurements and adjust themselves accordingly.

Of course, the figure size will by default not be appropriate for such an arrangement, and the content will clip:


\begin{examplefigure}{svg = true}
```julia
f = Figure()
for i in 1:5, j in 1:5
    Axis(f[i, j], width = 150, height = 150)
end
f
```
\end{examplefigure}


But like before we can call `resize_to_layout!` and the size will be corrected so no clipping occurs.


\begin{examplefigure}{svg = true}
```julia
set_theme!() # hide
resize_to_layout!(f)
f
```
\end{examplefigure}

# Scene tutorial

The scene constructor:

```julia
scene = Scene(;
    # clear everything behind scene
    clear = true,
    # the camera struct of the scene.
    visible = true,
    # ssao and light are explained in more detail in `Documetation/Lighting`
    ssao = Makie.SSAO(),
    # Creates lights from theme, which right now defaults to `
    # set_theme!(lightposition=:eyeposition, ambient=RGBf(0.5, 0.5, 0.5))`
    lights = Makie.automatic,
    backgroundcolor = :gray,
    resolution = (500, 500);
    # gets filled in with the currently set global theme
    theme_kw...
)
```

A scene is doing three things:

* holds a local theme, that gets applied to all plot objects in that scene
* manages the camera, projection and transformation matrices
* defines the window size. For sub-scenes, the child scene can have smaller window areas than the parent area.
* holds a reference to all window events

## Scenes and subwindows

With scenes, one can create subwindows. The window extends are given by a `Rect{2, Int}` and the position is always in window pixels and relative to the parent.

\begin{examplefigure}{}
```julia
using GLMakie, Makie
GLMakie.activate!()
scene = Scene(backgroundcolor=:gray)
subwindow = Scene(scene, px_area=Rect(100, 100, 200, 200), clear=true, backgroundcolor=:white)
scene
```
\end{examplefigure}

When using `Scenes` directly, one needs to manually set up the camera
and center the camera to the content of the scene
As described in more detail the camera section, we have multiple `cam***!` functions to set a certain projection and camera type for the scene.

\begin{examplefigure}{}
```julia
GLMakie.activate!() # hide
cam3d!(subwindow)
meshscatter!(subwindow, rand(Point3f, 10), color=:gray)
center!(subwindow)
scene
```
\end{examplefigure}

Instead of a white background, we can also stop clearing the background
to make the scene see-through, and give it an outline instead.
The easiest way to create an outline is, to make a sub scene with a projection that goes from 0..1 for the whole window.
To make a subscene with a certain projection type, Makie offers for each camera function a version without `!`, that will create a subscene, and apply the camera type.
We call the space that goes from 0..1 `relative` space, so `camrelative` will give this projection:

\begin{examplefigure}{}
```julia
GLMakie.activate!() # hide
subwindow.clear = false
relative_space = Makie.camrelative(subwindow)
# this draws a line at the scene window boundary
lines!(relative_space, Rect(0, 0, 1, 1))
scene
```
\end{examplefigure}

We can also now give the parent scene a more exciting background by using `campixel!` and plotting an image to the window:

\begin{examplefigure}{}
```julia
GLMakie.activate!() # hide
campixel!(scene)
w, h = size(scene) # get the size of the scene in pixels
# this draws a line at the scene window boundary
image!(scene, [sin(i/w) + cos(j/h) for i in 1:w, j in 1:h])
scene
```
\end{examplefigure}

We can fix this by translating the scene further back:

\begin{examplefigure}{}
```julia
GLMakie.activate!() # hide
translate!(scene.plots[1], 0, 0, -1000)
scene
```
\end{examplefigure}

We need a fairly high translation, since the far + near plane for `campixel!` goes from `-1000` to `1000`, while for `cam3d!` those get automatically adjusted to the camera parameters. Both end up in the same depthbuffer, transformed to the range `0..1` by the far & near plane, so to stay behind the 3d scene, it needs to be set to a high value.

With `clear = true` we wouldn't have this problem!

In GLMakie, we can actually take a look at the depthbuffer, to see how it looks now:

\begin{examplefigure}{}
```julia
GLMakie.activate!() # hide
screen = display(scene) # use display, to get a reference to the screen object
depth_color = GLMakie.depthbuffer(screen)
# Look at result:
f, ax, pl = heatmap(depth_color)
Colorbar(f[1, 2], pl)
f
```
\end{examplefigure}


## Window Events

Every scene also holds a reference to all global window events:
```julia:ex-scene
scene.events
```
\show{ex-scene}

We can use those events to e.g. move the subwindow. If you execute the below in GLMakie, you can move the sub-window around by pressing left mouse & ctrl:

```julia
on(scene.events.mouseposition) do mousepos
    if ispressed(subwindow, Mouse.left & Keyboard.left_control)
        subwindow.px_area[] = Rect(Int.(mousepos)..., 200, 200)
    end
end
```

## Projections and Camera

We've already talked a bit about cameras, but not really how it works.
Lets start from zero. By default, the scene x/y extends go from -1 to 1.
So, to draw a rectangle outlining the scene window, the following rectangle does the job:
\begin{examplefigure}{}
```julia
GLMakie.activate!() # hide
scene = Scene(backgroundcolor=:gray)
lines!(scene, Rect2f(-1, -1, 2, 2), linewidth=5, color=:black)
scene
```
\end{examplefigure}

this is, because the projection matrix and view matrix are the identity matrix by default, and Makie's unit space is what's called `Clip space` in the OpenGL world

```julia:ex-scene
cam = Makie.camera(scene) # this is how to access the scenes camera
```
\show{ex-scene}

One can change the mapping, to e.g. draw from -3 to 5 with an orthographic projection matrix:

\begin{examplefigure}{}
```julia
GLMakie.activate!() # hide
cam.projection[] = Makie.orthographicprojection(-3f0, 5f0, -3f0, 5f0, -100f0, 100f0)
scene
```
\end{examplefigure}

one can also change the camera to a perspective 3d projection:

\begin{examplefigure}{}
```julia
GLMakie.activate!() # hide
w, h = size(scene)
nearplane = 0.1f0
farplane = 100f0
aspect = Float32(w / h)
cam.projection[] = Makie.perspectiveprojection(45f0, aspect, nearplane, farplane)
# Now, we also need to change the view matrix
# to "put" the camera into some place.
eyeposition = Vec3f(10)
lookat = Vec3f(0)
upvector = Vec3f(0, 0, 1)
cam.view[] = Makie.lookat(eyeposition, lookat, upvector)
scene
```
\end{examplefigure}

## Interaction with Axis & Layouts

The Axis contains a scene, which has the projection set to make the coordinates go from `(x/y)limits_min ... (x/y)limits_max`. That's what we plot into.
Besides that, it's a normal scene, which we can use to create subscenes with smaller window size or a different projection.

So, we can use `camrelative` and friends to e.g. plot in the middle of the axis:

\begin{examplefigure}{}
```julia
GLMakie.activate!() # hide
figure, axis, plot_object = scatter(1:4)
relative_projection = Makie.camrelative(axis.scene);
scatter!(relative_projection, [Point2f(0.5)], color=:red)
# offset & text are in pixelspace
text!(relative_projection, "Hi", position=Point2f(0.5), offset=Vec2f(5))
lines!(relative_projection, Rect(0, 0, 1, 1), color=:blue, linewidth=3)
figure
```
\end{examplefigure}


## Transformations and Scene graph

So far we've been discussing only camera transformations of the scene.
In contrast, there are also scene transformations, or commonly referred to as world transformations.
To learn more about the different spaces, [learn opengl](https://learnopengl.com/Getting-started/Coordinate-Systems) offers some pretty nice explanations

The "world" transformation is implemented via the `Transformation` struct in Makie. Scenes and plots both contain these, so these types are considered as "Makie.Transformable".
The transformation of a scene will get inherited by all plots added to the scene.
An easy way to manipulate any `Transformable` is via these 3 functions:

{{doc translate!}}
{{doc rotate!}}
{{doc scale!}}

\begin{examplefigure}{}
```julia
GLMakie.activate!() # hide
scene = Scene()
cam3d!(scene)
sphere_plot = mesh!(scene, Sphere(Point3f(0), 0.5), color=:red)
scale!(scene, 0.5, 0.5, 0.5)
rotate!(scene, Vec3f(1, 0, 0), 0.5) # 0.5 rad around the y axis
scene
```
\end{examplefigure}

One can also transform the plot objects directly, which then adds the transformation from the plot object on top of the transformation from the scene.
One can add subscenes and interact with those dynamically.
Makie offers here what's usually referred to as a scene graph.

\begin{examplefigure}{}
```julia
GLMakie.activate!() # hide
translate!(sphere_plot, Vec3f(0, 0, 1))
scene
```
\end{examplefigure}

The scene graph can be used to create rigid transformations, like for a robot arm:

\begin{examplefigure}{}
```julia
GLMakie.activate!() # hide
parent = Scene()
cam3d!(parent)

# One can set the camera lookat and eyeposition, by getting the camera controls and using `update_cam!`
camc = cameracontrols(parent)
update_cam!(parent, camc, Vec3f(0, 8, 0), Vec3f(4.0, 0, 0))

s1 = Scene(parent, camera=parent.camera)
mesh!(s1, Rect3f(Vec3f(0, -0.1, -0.1), Vec3f(5, 0.2, 0.2)))
s2 = Scene(s1, camera=parent.camera)
mesh!(s2, Rect3f(Vec3f(0, -0.1, -0.1), Vec3f(5, 0.2, 0.2)), color=:red)
translate!(s2, 5, 0, 0)
s3 = Scene(s2, camera=parent.camera)
mesh!(s3, Rect3f(Vec3f(-0.2), Vec3f(0.4)), color=:blue)
translate!(s3, 5, 0, 0)
parent
```
\end{examplefigure}

\begin{examplefigure}{}
```julia
GLMakie.activate!() # hide
# Now, rotate the "joints"
rotate!(s2, Vec3f(0, 1, 0), 0.5)
rotate!(s3, Vec3f(1, 0, 0), 0.5)
parent
```
\end{examplefigure}

With this basic principle, we can even bring robots to life :)
[Kevin Moerman](https://github.com/Kevin-Mattheus-Moerman) was so nice to supply a Lego mesh, which we're going to animate!
When the scene graph is really just about a transformation Graph, one can use the Transformation struct directly, which is what we're going to do here.
This is more efficient and easier than creating a scene for each model.
Let's use WGLMakie with it's offline export feature, to create a plot with sliders to move the parts, that keeps working in the browser:

\begin{showhtml}{}
```julia
using WGLMakie, JSServe
WGLMakie.activate!()
Page(offline=true, exportable=true)
```
\end{showhtml}

\begin{showhtml}{}
```julia
using MeshIO, FileIO, GeometryBasics

colors = Dict(
    "eyes" => "#000",
    "belt" => "#000059",
    "arm" => "#009925",
    "leg" => "#3369E8",
    "torso" => "#D50F25",
    "head" => "yellow",
    "hand" => "yellow"
)

origins = Dict(
    "arm_right" => Point3f(0.1427, -6.2127, 5.7342),
    "arm_left" => Point3f(0.1427, 6.2127, 5.7342),
    "leg_right" => Point3f(0, -1, -8.2),
    "leg_left" => Point3f(0, 1, -8.2),
)

rotation_axes = Dict(
    "arm_right" => Vec3f(0.0000, -0.9828, 0.1848),
    "arm_left" => Vec3f(0.0000, 0.9828, 0.1848),
    "leg_right" => Vec3f(0, -1, 0),
    "leg_left" => Vec3f(0, 1, 0),
)

function plot_part!(scene, parent, name::String)
    # load the model file
    m = load(assetpath("lego_figure_" * name * ".stl"))
    # look up color
    color = colors[split(name, "_")[1]]
    # Create a child transformation from the parent
    child = Transformation(parent)
    # get the transformation of the parent
    ptrans = Makie.transformation(parent)
    # get the origin if available
    origin = get(origins, name, nothing)
    # center the mesh to its origin, if we have one
    if !isnothing(origin)
        centered = m.position .- origin
        m = GeometryBasics.Mesh(meta(centered; normals=m.normals), faces(m))
        translate!(child, origin)
    else
        # if we don't have an origin, we need to correct for the parents translation
        translate!(child, -ptrans.translation[])
    end
    # plot the part with transformation & color
    return mesh!(scene, m; color=color, transformation=child)
end

function plot_lego_figure(s, floor=true)
    # Plot hierarchical mesh and put all parts into a dictionary
    figure = Dict()
    figure["torso"] = plot_part!(s, s, "torso")
        figure["head"] = plot_part!(s, figure["torso"], "head")
            figure["eyes_mouth"] = plot_part!(s, figure["head"], "eyes_mouth")
        figure["arm_right"] = plot_part!(s, figure["torso"], "arm_right")
            figure["hand_right"] = plot_part!(s, figure["arm_right"], "hand_right")
        figure["arm_left"] = plot_part!(s, figure["torso"], "arm_left")
            figure["hand_left"] = plot_part!(s, figure["arm_left"], "hand_left")
        figure["belt"] = plot_part!(s, figure["torso"], "belt")
            figure["leg_right"] = plot_part!(s, figure["belt"], "leg_right")
            figure["leg_left"] = plot_part!(s, figure["belt"], "leg_left")

    # lift the little guy up
    translate!(figure["torso"], 0, 0, 20)
    # add some floor
    floor && mesh!(s, Rect3f(Vec3f(-400, -400, -2), Vec3f(800, 800, 2)), color=:white)
    return figure
end
App() do session
    s = Scene(resolution=(500, 500))
    cam3d!(s)
    figure = plot_lego_figure(s, false)
    bodies = [
        "arm_left", "arm_right",
        "leg_left", "leg_right"]
    sliders = map(bodies) do name
        slider = if occursin("arm", name)
            JSServe.Slider(-60:4:60)
        else
            JSServe.Slider(-30:4:30)
        end
        rotvec = rotation_axes[name]
        bodymesh = figure[name]
        on(slider) do val
            rotate!(bodymesh, rotvec, deg2rad(val))
        end
        DOM.div(name, slider)
    end
    center!(s)
    JSServe.record_states(session, DOM.div(sliders..., s))
end
```
\end{showhtml}

Finally, lets let him walk and record it as a video with the new, experimental ray tracing backend.

Note: RPRMakie is still not very stable and rendering out the video is quite slow on CI, so the shown video is prerendered!

```julia
using RPRMakie
# iterate rendering 200 times, to get less noise and more light
RPRMakie.activate!(iterations=200)

radiance = 50000
# Note, that only RPRMakie supports `EnvironmentLight` so far
lights = [
    EnvironmentLight(1.5, rotl90(load(assetpath("sunflowers_1k.hdr"))')),
    PointLight(Vec3f(50, 0, 200), RGBf(radiance, radiance, radiance*1.1)),
]
s = Scene(resolution=(500, 500), lights=lights)
cam3d!(s)
c = cameracontrols(s)
c.near[] = 5
c.far[] = 1000
update_cam!(s, c, Vec3f(100, 30, 80), Vec3f(0, 0, -10))
figure = plot_lego_figure(s)

rot_joints_by = 0.25*pi
total_translation = 50
animation_strides = 10

a1 = LinRange(0, rot_joints_by, animation_strides)
angles = [a1; reverse(a1[1:end-1]); -a1[2:end]; reverse(-a1[1:end-1]);]
nsteps = length(angles); #Number of animation steps
translations = LinRange(0, total_translation, nsteps)

Makie.record(s, "lego_walk.mp4", zip(translations, angles)) do (translation, angle)
    #Rotate right arm+hand
    for name in ["arm_left", "arm_right",
                            "leg_left", "leg_right"]
        rotate!(figure[name], rotation_axes[name], angle)
    end
    translate!(figure["torso"], translation, 0, 20)
end
```
~~~
<video autoplay controls src="/assets/lego_walk.mp4">
</video>
~~~
# MakieRecipes

[![build status](https://github.com/JuliaPlots/MakieRecipes.jl/workflows/CI/badge.svg)](https://github.com/JuliaPlots/MakieRecipes.jl/actions)
[![Docs][docs-img]][docs-url]

[docs-img]: https://img.shields.io/badge/docs-dev-blue.svg
[docs-url]: http://juliaplots.org/MakieRecipes.jl/dev/

Extending Makie to support Plots.jl recipes!
