# NoisySignalIntegration.jl

*A tool to determine uncertainty in numeric integrals of noisy x-y data.*

`NoisySignalIntegration` implements a method to determine the uncertainty in
numeric integrals of noisy x-y data on the basis of a Monte-Carlo process.  It
can include uncertainty due to noise, baseline subtraction, and placement in
integration bounds.  To do this, the integration is repeated many times while
the noise of the data, baseline, and integration bounds are varied based on a
noise model and user supplied probability distributions.

To view the documentation, click the badge below:

[![Documentation, latest](https://img.shields.io/badge/docs-latest-blue.svg)](https://nluetts.github.io/NoisySignalIntegration.jl/dev/)

## Installation

The package is not yet registered in Julia's general package registry, you have to install it directly from this Github repository.

To install it for your project, enter the package mode in the Julia REPL (press `]`) and type:

```
add https://github.com/nluetts/NoisySignalIntegration.jl
```

While still in package mode, you can type

```
test NoisySignalIntegration
```

to run the package's unit tests.

## Getting Started

Check out the [documentation](https://nluetts.github.io/NoisySignalIntegration.jl/dev/) to learn how to use the package.

If you don't have a local Julia installation, you can test the package on [mybinder.org](https://mybinder.org).
Click the badge below and open one of the example notebooks:

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/nluetts/NSI-Binder/HEAD?filepath=example-1.ipynb&urlpath=lab)

## Contributing and Support

If you have problems with the package, submit an [issue](https://github.com/nluetts/NoisySignalIntegration.jl/issues).
Feel free to fork the project and open a pull request if you would like to contribute.# NoisySignalIntegration.jl

*NoisySignalIntegration.jl -- A tool to determine uncertainty in numeric
integrals of noisy x-y data.*

`NoisySignalIntegration` implements a method to determine the uncertainty in
numeric integrals of noisy x-y data on the basis of a Monte-Carlo process.  It
can include uncertainty due to noise, baseline subtraction, and placement in
integration bounds.  To do this, the integration is repeated many times while
the noise of the data, baseline, and integration bounds are varied based on a
noise model and user supplied probability distributions.

A predecessor of this package was originally intended to estimate uncertainties
of band signals in FTIR spectra (see [G. Karir et al.,
2019](https://doi.org/10.1039/C9CP00435A)), which is reflected in the example
given in the [Usage Guide](@ref).

**Table of Contents**

```@contents
Pages = ["overview.md", "guide.md", "examples.md", "baseline.md", "API.md", "internals.md"]
Depth = 4
```
NoisySignalIntegration.jl is available from [Github.com](https://github.com/nluetts/NoisySignalIntegration.jl).

To install it for your project, enter the package mode in the Julia REPL (press `]`) and type:

```
add https://github.com/nluetts/NoisySignalIntegration.jl
```

While still in package mode, you can type

```
test NoisySignalIntegration
```

to run the package's unit tests.

Note that the package requires Julia v1.5 or above.# Usage Guide

As a more detailed usage example, we will go through the analysis of a simulated
FTIR spectrum.

!!! info "A note on plotting"
    `NoisySignalIntegration` provides several
    "[recipes](http://docs.juliaplots.org/latest/recipes/)" for
    [`Plots.jl`](http://docs.juliaplots.org/latest/) to easily plot the various
    (interim) results. Often, merely calling `plot()` and passing in data types
    from `NoisySignalIntegration` will work. Examples are included in this
    guide.

Suppose our spectrum looks like the following simulation:

```@example FTIR
using Distributions, NoisySignalIntegration, Plots
using Random: seed!

spectrum = NoisySignalIntegration.testdata_1()
plot(spectrum, label="simulated spectrum")
```

In order two apply the `NoisySignalIntegration` uncertainty analysis, we must perform 4 basic steps:

1. From the spectrum, crop the region that contains the signals and the region
   that contains a representative sample of the noise
1. Characterize the noise (to be able to simulate it in the Monte-Carlo draws) 
1. Set integration bounds and their associated uncertainties
1. Run the [`mc_integrate`](@ref) function

## [Cropping the spectrum](@id crop)

Let's start by dividing the spectrum into the bands we want to integrate and the
noise we want to analyse. We can do this by using the [`crop`](@ref) function.

```@example FTIR
slice_bands = crop(spectrum,  5.0,  40.0)
slice_noise = crop(spectrum, 40.0, 100.0)

plot(slice_bands; label="bands")
plot!(slice_noise; label="noise")
```

## Noise analysis

The spectrum has a quite considerable baseline which constitutes a problem when
analysing the noise. To prepare the noise spectrum `slice_noise` for analysis,
we create a [`NoiseSample`](@ref) object. Upon construction of the `NoiseSample`
object, a polynomial is fitted and subtracted to remove the baseline:

```@example FTIR
noise = NoiseSample(slice_noise, 3) # 3 = remove third order polynomial baseline
nothing # hide
```

Plotting the original slice and the `NoiseSample` object shows the data after
baseline removal:

```@example FTIR
plot(slice_noise, label="cropped noise")
plot!(noise, label="NoiseSample object")
```

In order to simulate the noise, we must determine its characteristics. A model
is retrieved by fitting the estimated autocovariance:

```@example FTIR
nm = fit_noise(noise)

# plot the fitting result:
plotautocovfit(noise, nm);
# create zoomed inset:
lens!([0, 1.5], [-1e-3, 3e-3], inset = (1, bbox(0.3, 0.3, 0.3, 0.3)))
```

Plotting the model next to the noise object is an important sanity check to
verify that the fitting yielded a sensible estimate and that generated noise
samples do mimic the experimental noise.

```@example FTIR
plot(noise, nm)
```

!!! info "plotting more samples"
    They keyword `draws` controls how many random generated noise draws are
    plotted:
    
    ```julia
    plot(noise, nm; draws=5) # draw 5 instead of 3 (default) noise draws
    ```

    This also works when you [plot Monte-Carlo draws](@ref plot_mc_draws).

## Preparing the spectrum for integration

Now that we have a noise model, we can generate an [`UncertainCurve`](@ref).
An `UncertainCurve` holds random draws of the original spectrum plus noise:

```@example FTIR
uncertain_spectrum = add_noise(slice_bands, nm, 50_000) # generate 50_000 random samples
nothing # hide
```

If we plot the `uncertain_spectrum`, we get a ribbon plot showing a 95%
confidence band:

```@example FTIR
plot(uncertain_spectrum)
```

We can also plot single draws by using the `mcplot()` function from
`MonteCarloMeasurements.jl`:

```@example FTIR
using MonteCarloMeasurements

mcplot(uncertain_spectrum; draws=20)
```


## [Integration bounds](@id bounds_guide)

`NoisySignalIntegration` deals with uncertainty in placing integration bounds by
expressing each bound by one ore more probability distributions. Any continuous,
univariate distribution from
[Distributions.jl](https://github.com/JuliaStats/Distributions.jl) can be used
to define integration bounds.

To define integration bounds, the `UncertainBound` object is used. There are
several options available to create an `UncertainBound`:

1. Passing two distributions
1. Passing a position, a distribution, and an `UncertainCurve`
1. Passing a vector of positions, a distribution, and an `UncertainCurve`

### Defining an `UncertainBound` using a start and end point

If two distributions are passed to `UncertainBound()`, they will be interpreted
as start and end points for the integration window with the uncertainty of these
points being expressed by the spread of the distributions.

For example, let's say we want to integrate the right peak of our simulated
spectrum:

```@example FTIR
plot(crop(spectrum, 20, 40), label="right peak")
plot!([27, 32], [1.3, 1.3]; markershape=:cross, label="integration interval")
```

It looks like integrating from about 27 to 32 would be appropriate, but there is
some doubt of the exact location of the integration bounds. Perhaps a reasonable
estimate is that the left bound falls in the range from 26 to 27 and the right
bound in the range from 32 to 33. This would be expressed with a
`UncertainBound` that is defined using two uniform distributions:

```@example FTIR
lrb = UncertainBound(Uniform(26, 27), Uniform(32, 33)) # 10 k samples by default
nothing # hide
```

Upon creation of the `UncertainBound` object, pseudo random samples of the
integration start and end point are drawn. If we do not provide the number of
samples, it will default to 10 000. We can plot the bound as a histogram to see
the distribution of the start and end point:

```@example FTIR
histogram(lrb; label=["start" "end"])
```

The uniform distribution is of course a bit of an awkward choice, because its
probability density suddenly drops to 0, which perhaps does not model one's
belief about the position of the integration bounds very well.

Due to the central limit theorem and the general applicability of the normal
distribution, it is often a natural choice when dealing with uncertainties:

```@example FTIR
lrb_normal = UncertainBound(Normal(26.5, 0.5), Normal(32.5, 0.5), 12_000) # we draw 12_000 samples, just to illustrate how it works

histogram(lrb_normal; label=["start" "end"])
```

However, in this particular case of describing the uncertainty of integration
bounds, the tails of the normal distribution are problematic, because they lead
to occasional extreme values of the integration interval, which would not seem
realistic.

A compromise between the uniform and normal distribution is a scaled and shifted
beta(2, 2) distribution. Its shape resembles the shape of the normal
distribution but it is missing the tails. Since a scaled and shifted beta
distribution does not ship with
[Distributions.jl](https://github.com/JuliaStats/Distributions.jl),
`NoisySignalIntegration` includes the function [`scale_shift_beta`](@ref)`(α, β,
a, b)` which can be used to generate a beta(α, β) distribution that has a
support in the interval `a` to `b`.

Again, a demonstration may help to explain (we keep the normal distribution for
the right bound so we can compare the distributions):

```@example FTIR
lrb_beta_normal = UncertainBound(scale_shift_beta(2, 2, 26, 27), Normal(32.5, 0.5))

histogram(lrb_beta_normal; label=["start" "end"], normalize=true)
```

### Defining an `UncertainBound` using a position, width, and `UncertainCurve` (symmetric bands)

For some spectra, one can assume that a band is more or less symmetric. In this
case, it may be better to define an integration window not by a start and end
point but merely by a width, and integrate the band symmetrically around its
peak position ± half this width.

To accomplish this, one has to construct an `UncertainBound` object by passing a
`position` (the peak position of the symmetic band), a distribution that
describes the uncertainty in the integration window width (here
`width_distribution`), and an `UncertainCurve` that holds samples of the
spectrum (here `uncertain_spectrum`):

```@example FTIR
position = 15.0
 # widths will fall in the range 2 to 3, with a maximum at 2.5
width_distribution = scale_shift_beta(2, 2, 2, 3)
# define a "width bound"
wb = UncertainBound(position, width_distribution, uncertain_spectrum)
nothing # hide
```

From the provided data, the `UncertainBound` object is created as follows:
- The width distribution is used to draw as many random samples of the
  integration window width $w$ as the `uncertain_spectrum` contains spectral
  samples
- For each spectral sample in `uncertain_spectrum`, the peak position $p_x$ in
  the range `position` ± $\frac{w}{2}$ is retrieved
- The peak position $p_x$ is used to define the start and end point of the
  integration window for each spectral sample, $p_x - \frac{w}{2}$ and $p_x +
  \frac{w}{2}$

Therefore, what is stored in the `UncertainBound` object are again start and end
points for the integration. We can verify this by plotting another histogram:

```@example FTIR
histogram(wb; label=["start" "end"], normalize=true, linewidth=0)
```

The crucial difference compared to the bound defined from two distributions is
that the start and end points are now placed symmetrically with respect to the
band's peak position. The `UncertainBound` now "follows" the peak in each
Monte-Carlo draw, so to speak.

### Defining an `UncertainBound` using several positions, a width, and `UncertainCurve` (several symmetric bands with same width)

If, for example from a physical argument, we can say that two bands should have
the same width, we can constrain our `UncertainBound`s even further: we can
create several bounds that share the exact same integration window width in each
draw.

All we have to do is to provide the constructor of `UncertainBound` not with a
single position, but with an array of several positions:

```@example FTIR
positions = [15.0, 30.0]

wb_1, wb_2 = UncertainBound(positions, width_distribution, uncertain_spectrum)
nothing # hide
```

Note that the constructor will then return an array of `UncertainBound` objects
which we unpacked into the variables `wb_1` and `wb_2` in the example above.

The histograms of the start and end points looks like this:

```@example FTIR
histogram( wb_1; label=["start 1" "end 1"], normalize=true, linewidth=0)
histogram!(wb_2; label=["start 2" "end 2"], normalize=true, linewidth=0)
```

It is not obvious from the histograms that the widths of the integration windows
stored in `wb_1` and `wb_2` are identical, so we calculate and print them here
manually to prove this:

```@example FTIR
for i in 1:10
    l1 = wb_1.left.particles[i]
    l2 = wb_2.left.particles[i]
    r1 = wb_1.right.particles[i]
    r2 = wb_2.right.particles[i]

    println("draw $i: width 1 = ", r1 - l1, " width 2 = ", r2 - l2)
end
```

!!! warning "Watch out for the support of your width distribution"
    Note that the distribution that you pass to `UncertainBound` along with a
    position/positions must not allow for negative values (i.e. its support must
    end before 0). Keep in mind that a normal distribution, for
    example, has support from -∞ to ∞, so it is a poor choice here.

## [Plotting Monte-Carlo draws](@id plot_mc_draws)

To verify that the integration windows and derived integrals are sensible, it is
a good idea to plot a few draws and integrals before running the full
Monte-Carlo algorithm. We can do so by passing an `UncertainCurve` and an array
of `UncertainBound`s to the plot function:

```@example FTIR
plot(uncertain_spectrum, [wb_1, wb_2]; size=(400, 500), local_baseline=true)
```

We can see from the plot that our estimate for the width of the peaks was
perhaps a bit too small, so we retry:

```@example FTIR
width_distribution = scale_shift_beta(2, 2, 3, 4) # width will fall in the range [3, 4]
wb_1, wb_2 = UncertainBound(positions, width_distribution, uncertain_spectrum)

plot(uncertain_spectrum, [wb_1, wb_2]; size=(400, 500), local_baseline=true)
```

An alternative to the plot function is the [`animate_draws`](@ref) function which allows you
to visualize the Monte-Carlo draws in a gif animation:

```@example FTIR
NoisySignalIntegration.animate_draws(
    uncertain_spectrum, [wb_1, wb_2];
    size=(300, 150),
    local_baseline=true
)
```

## Running the integration algorithm

The integration is performed with the function [`mc_integrate`](@ref). We have
to pass in the uncertain spectrum and integration bounds. Since we pass in two
integration bounds, we retrieve two areas:

```@example FTIR
area_1, area_2 = mc_integrate(uncertain_spectrum, [wb_1, wb_2]; local_baseline=true)
```

We can look at the histogram of the integrals:

```@example FTIR
histogram([area_1.particles, area_2.particles]; label=["band area 1" "band area 2"])
```

Or of the peak area ratio, simply by calculating with the retrieved areas:

```@example FTIR
ratio = area_1 / area_2
histogram(ratio.particles; label="band area ratio (band 1/band 2)")
```

We see that the histogram of peak area ratio peaks around 0.5, which is what we
put into the simulation of the spectrum.

We can use some basic statistical functions to characterize the result:

```@example FTIR
using StatsBase: mean, std, percentile

mean(ratio)
```

```@example FTIR
std(ratio)
```

```@example FTIR
percentile(ratio, 2.5)
```

```@example FTIR
percentile(ratio, 50)
```

```@example FTIR
percentile(ratio, 97.5)
```

We find that, considering the noise and the uncertainty in the integration
bounds, we end up with a 95% uncertainty interval of roughly 0.4 to 0.7.

!!! info "Sensitivity analysis"
    It is perfectly valid to create several `UncertainBound`s for *one and the
    same band* and feed them into `mc_integrate()`, e.g. to perform a sensitivity
    analysis on how much the result depends on the kind and parameters of the
    bounds.

You find more usage examples on the next page. In particular, check out the
[error propagation example](@ref example_propagation) to see how to proceed with
uncertainty calculations with the retrieved areas using
[`MonteCarloMeasurements.jl`](https://github.com/baggepinnen/MonteCarloMeasurements.jl).# API Reference

```@contents
Pages = ["API.md"]
Depth = 4
```

## Types

```@docs
Curve

NoiseSample

UncertainCurve

GaussianNoiseModel

MvGaussianNoiseModel

UncertainBound
```

## Functions

### Manipulation of Curves

```@docs
crop

stitch

add_noise
```

### Noise analysis

```@docs
fit_noise

plotautocovfit

```

### Statistics

```@docs
scale_shift_beta
```

### Integration

```@docs
mc_integrate

trapz
```

## Macros

```@docs
@samples
```

## Misc

```@docs
animate_draws
```# Baseline Handling

There are several ways to handle baseline correction when working with the package.
The easiest method is to use the build-in local baseline correction which
assumes a linear baseline between the start and end point of the integration
window. It is envoked by using the keyword argument `subtract_baseline`
(deprecated in v0.2) or `local_baseline`. The difference of these methods is discussed
below.

Otherwise, one can subtract a baseline from the data in a preprocessing step.
This can be done either before or after generating an [`UncertainCurve`](@ref).
If a baseline is subtracted from the `UncertainCurve`, it is possible to
account for uncertainty in the baseline correction, e.g. by subtracting
baselines generated using a Gaussian process.

## Build-in

Local linear baseline subtraction can be achieved by passing
`local_baseline=true` to the [`mc_integrate`](@ref) function. To visualize the
integrated area, the same keyword argument can be passed to the `plot()`
function when [plotting Monte-Carlo draws](@ref plot_mc_draws).

The keyword argument `subtract_baseline=true` is also supported, but its use is
deprecated. The difference of `local_baseline` and `subtract_baseline` can be
visualized when animating draws of curves and integration bound samples (using
`Plots.@animate`): 

```@eval
using NoisySignalIntegration
using Plots
using Random: seed!

seed!(1)
nsi = NoisySignalIntegration

let
    n = 20
    c = crop(nsi.testdata_1(), 0, 50)
    uc = add_noise(c, GaussianNoiseModel(0.1))
    ubleft = UncertainBound(15.0, scale_shift_beta(2.0, 2.0, 4.0, 5.0), uc)
    ubright = UncertainBound(30.0, scale_shift_beta(2.0, 2.0, 6.0, 7.0), uc)

    spany = [
        mm(curve.y for curve ∈ [nsi.get_draw(i, uc) for i ∈ 1:n]) |> mm
        for mm in (minimum, maximum)
    ]
    spany = (spany[1]*0.9, spany[2]*1.1)
    anim = @animate for i in 1:n
        kw = Dict(:ylim => spany, :legend => :topleft)
        p1 = plot(nsi.get_draw(i, uc), [ubleft, ubright], i; subtract_baseline=true, label=["subtract_baseline" "" ""], kw...)
        p2 = plot(nsi.get_draw(i, uc), [ubleft, ubright], i; local_baseline=true, label=["local_baseline" "" ""], kw...)
        plot(p1, p2; layout=(2, 1), size=(400, 300))
    end
    gif(anim, "baseline_anim.gif", fps=5)
    nothing
end
```

![build-in baseline handling animation](baseline_anim.gif)

As the figure above shows, the baseline varies considerably more when
`subtract_baseline` is used, compared to `local_baseline`. The former uses the
exact start and end points of the integration window in a particular draw while
the latter uses the start and end point distributions to determine a weighted
average for y-values at integration bounds of a particular draw. Especially at high noise
levels, rather extreme local baseline estimates can follow from using
`subtract_baseline` and the overall uncertainty may be overestimated. It is
thus recommended to use the `local_baseline` method.

## Preprocessing

### Simple baseline subtraction

The data can simply be preprocessed before a [`Curve`](@ref) object is created.
For example, one can mask the signals that shall be integrated (using a filter
or calls to the [`crop`](@ref) and [`stitch`](@ref) functions) and fit a
[polynomial](https://github.com/JuliaMath/Polynomials.jl) or
[smoothing spline](https://github.com/nignatiadis/SmoothingSplines.jl) to be
subtracted:



```@setup simple_baseline
using NoisySignalIntegration, Plots
crv = NoisySignalIntegration.testdata_1()
```

```@example simple_baseline
using Polynomials, SmoothingSplines

polyfit, splinefit = let
    no_signals = stitch(crop(crv, 0, 12), crop(crv, 32, 300)) # remove signals
    x = no_signals.x
    y = no_signals.y
    pfit = Polynomials.fit(x, y, 5)
    sfit = SmoothingSplines.fit(SmoothingSpline, x, y, 2000.0)
    pfit, sfit
end

crv_baseline_corrected = crv - predict(splinefit, crv.x)

plot(crv; label="data", alpha=0.3, color=:black, legend=:outertopright, size=(800, 400))
plot!(polyfit, 0, 100; label="polynomial baseline fit")
plot!(crv.x, predict(splinefit, crv.x); label="spline baseline fit")
plot!(crv_baseline_corrected; label="data - spline baseline")

```

### Uncertain baseline

It may be desirable to account for uncertainty in the subtracted baseline.
One way to achieve this is to fit a [Gaussian process](https://github.com/STOR-i/GaussianProcesses.jl) to the dataset while
excluding signals. For example, consider a dataset with a relatively broad
band:

```@setup ubaseline
using NoisySignalIntegration, Plots, Polynomials

baseline = Polynomial([1, -0.003, 1e-5, 3e-7]).(-150:0.5:150)
nm = MvGaussianNoiseModel(0.5, 0.1, 1.5);

crv = NoisySignalIntegration.generate_testdata(
    0:0.5:300,
    [(15, 40, 5), (30, 140, 10), (10, 150, 15), (5, 155, 8), (25, 170, 10)],
    nm;
    seedvalue=42,
    baseline=baseline
)
```

```@example ubaseline
plot(crv, label=nothing)
```

A Gaussian process can be used to approximate the baseline, while allowing for
higher uncertainty in the regions where the course of the baseline is masked by
signals:


```@example ubaseline
using GaussianProcesses, MonteCarloMeasurements

gp = let
    # cut out signals
    no_signals = stitch(crop(crv, 0, 25), crop(crv, 55, 110), crop(crv, 200, 300))
    # fit Gaussian process
    GP(no_signals.x, no_signals.y, MeanZero(), SE(4.0, 0.1))
end

# draw random baselines ...
fs = rand(gp, crv.x, 250)

# ... and plot
plot(crv.x, fs; alpha=0.05, color=:red, dashed=:dashed, legend=nothing)
plot!(crv)
```

The "amplitude" of the baseline's uncertainty in the range of the signals
can be tweaked by the hyperparameters of the Gaussian process, in particular
the autocorrelation length scale of the covariance function (in the example
`SE()`, the "squared exponential" function, i.e. Gaussian covariance function,
with a length scale of 4.0 units).


The baseline samples can be subtracted from an `UncertainCurve` (here `ucrv`) as follows:

```julia
# make sure the number of samples matches, here 100 000 samples
ubaseline = rand(gp, crv.x, 100_000) |> transpose |> collect |> Particles
ucrv_baseline_corrected = UncertainCurve(ucrv.x, ucrv.y - ubaseline)
```

The resulting `UncertainCurve` now includes not only uncertainty due to noise
but also due to the baseline correction.
# Internals

Upon creation, the data types [`UncertainCurve`](@ref) and [`UncertainBound`](@ref) use the package [MonteCarloMeasurements.jl](https://github.com/baggepinnen/MonteCarloMeasurements.jl) to generate random samples. Furthermore, integration results are returned as [Particle](https://baggepinnen.github.io/MonteCarloMeasurements.jl/stable/api/#MonteCarloMeasurements.Particles) objects from MonteCarloMeasurements which enables you to use MonteCarloMeasurements for further [error propagation](@ref example_propagation). Check out the [MonteCarloMeasurements.jl documentation](https://baggepinnen.github.io/MonteCarloMeasurements.jl/stable/) to learn more about the sampling process and uncertainty calculations.# Case studies

## Raman spectra

Raman spectra are typically measured with a CCD camera where there is no obvious autocorrelation of the noise.
The integration works in general as outlined in the [Usage Guide](@ref). The main difference is the analysis 
and generation of noise. Instead of a multivariate `MvGaussianNoiseModel`, we have to use a `GaussianNoiseModel`:

```@example Raman
# simulate a Raman spectrum with two bands

using Plots: plot, plot!, histogram
using Statistics
using MonteCarloMeasurements

using NoisySignalIntegration


spectrum = NoisySignalIntegration.testdata_2()
plot(spectrum, label="simulated spectrum")
```

```@example Raman
# crop the bands and the noise sample

bands = crop(spectrum, 10, 40) - 1 # we can directly subtract 1 from the curve to remove the offset
noise = NoiseSample(crop(spectrum, 40, 100), 3)

plot(bands; label="bands")
plot!(noise; label="noise sample")
```

```@example Raman
# build the noise model

nm = noise |> std |> GaussianNoiseModel

uncertain_spectrum = add_noise(bands, nm)

plot(uncertain_spectrum)
```

```@example Raman
# plot some MC draws to check

bounds = UncertainBound([15., 30.], scale_shift_beta(2, 2, 4.5, 5), uncertain_spectrum)

plot(uncertain_spectrum, bounds; size=(500, 600))
```

```@example Raman
# integrate the bands and calculate the ratio

area1, area2 = mc_integrate(uncertain_spectrum, bounds)

histogram(area1/area2)
```

## Mexican hat wavelet

This example is more of a basic sanity check for the integration algorithm.
The integration of a mexican hat wavelet should yield an integral of 0.
For a noisy curve with uncertain integration bounds, the results should still
be compatible with an integral of 0, i.e. the 95% confidence interval of the
resulting area should include the value 0. If we decrease the noise and the
uncertainty in the integration bounds, the result should come closer and closer
to 0 (the area distribution should become narrower while still covering the 
value 0).

We go once through the regular workflow. First, we generate the test data:

```@example mh
using NoisySignalIntegration
using Plots
using Random: seed!

seed!(7) # seed for reproducibility

function mexican_hat_curve(noise_model)
    x = -20:0.1:100
    y = @. (1 - x^2)*exp(-x^2/2)
    mh = Curve(x, y)
    uncertain_mh = add_noise(mh, noise_model, 1)
    return NoisySignalIntegration.get_draw(1, uncertain_mh)
end

noisy_mh_curve = mexican_hat_curve(MvGaussianNoiseModel(0.1, 0.05, 0.5))

plot(noisy_mh_curve; label="noisy mexican hat")
```

Analysis of noise:

```@example mh
noise = NoiseSample(crop(noisy_mh_curve, 5, 100))
plot(noise; label="noise sample")
```

```@example mh
noise_model = fit_noise(noise)
```

```@example mh
plotautocovfit(noise, noise_model)
```

Check that generated noise samples look realistic:

```@example mh
plot(noise, noise_model; size=(500, 600))
```

Definition of integration bounds:

```@example mh
uncertain_mh_curve = add_noise(noisy_mh_curve, noise_model)
# definition using distributions for start and end point
bnd = UncertainBound(scale_shift_beta(2, 2, -5, -4), scale_shift_beta(2, 2, 4, 5))
plot(uncertain_mh_curve, bnd; size=(500, 600), xlim=(-25, 25))
```

Integration:

```@example mh
area = mc_integrate(uncertain_mh_curve, bnd)

histogram(area; label="area")
```

As expected, the mean of the derived areas is close to zero.

Now we look what happens if we successively decrease the noise level and
uncertainty in the integration window position. We define a function
that covers the complete integration workflow for a certain scaling factor
`f` that controls the noise level and variability of integration windows:

```@example mh
function analyze_noisy_mh(f)
    noisy_mh_curve = mexican_hat_curve(MvGaussianNoiseModel(0.01, 0.05*f, 0.1))
    noise = NoiseSample(crop(noisy_mh_curve, 5, 100))
    noise_model = fit_noise(noise)
    println(noise_model)
    uncertain_mh_curve = add_noise(noisy_mh_curve, noise_model)
    bnd = UncertainBound(scale_shift_beta(2, 2, -4.5-0.5*f, -4.5+0.5*f), scale_shift_beta(2, 2, 4.5-0.5*f, 4.5+0.5*f))
    return mc_integrate(uncertain_mh_curve, bnd)
end
```

We apply the analysis function to several `f` factors and plot the histograms
of the resulting area distributions:

```@example mh
seed!(2)

areas = [analyze_noisy_mh(f) for f in (1, 0.5, 0.25, 0.05)]

histogram(
    [a.particles for a in areas];
    label=["f = 1" "f = 0.5" "f = 0.25" "f = 0.05"],
    normalize=true,
    alpha=0.3,
    linetype=:stephist,
    fill=true,
    xlim=(-1, 1)
)
```

As we can see, with decreasing uncertainty, the result comes closer and closer to 0.

## [Error propagation](@id example_propagation)

Since integration results are returned as [Particle](https://baggepinnen.github.io/MonteCarloMeasurements.jl/stable/api/#MonteCarloMeasurements.Particles) objects from the [MonteCarloMeasurements.jl](https://github.com/baggepinnen/MonteCarloMeasurements.jl) package (see [Internals](@ref)), calculating combined uncertainties, for example when deriving abundance ratios from band integrals and calculated intensities, is rather simple. Consider the following spectrum with four bands:

```@example propagation
# simulate an FTIR spectrum with 4 bands

using Plots
using NoisySignalIntegration

spectrum = NoisySignalIntegration.testdata_4()
plot(spectrum, label="simulated spectrum")
```

Let's assume that these bands correspond to 4 different chemical species
(species A, B, C and D from left to right) and we know from quantum chemical
calculations what the intensity of each band should be per unit of substance.
With this information, we can calculate abundance ratios $R_{xy}$ likes so:

```math
R_{xy} = \frac{A_x / I_x}{A_y / I_y}
```

Where $A$ is the band integral and $I$ the calculated intensity for species $x$
or $y$, respectively.

The uncertainty calculation for the abundance ratio is not straightforward if
the distribution of the integrals is asymmetric (which can easily be the case
for larger uncertainties). Using `MonteCarloMeasurement` makes it rather simple,
however.

First, we integrate with `NoisySignalIntegration`:

```@example propagation
# crop spectrum and noise
bands = crop(spectrum, 0, 100)
noise = NoiseSample(crop(spectrum, 100, 200), 3)
# retrieve noise model
nm = fit_noise(noise)
# prepare spectral samples
uspec = add_noise(bands, nm)
# declare integration bounds (several symmetric bands with same width)
bnds = UncertainBound([15., 30., 60., 85.], scale_shift_beta(2, 2, 4, 6), uspec)
# integrate
areas = mc_integrate(uspec, bnds; local_baseline=true)
A_A, A_B, A_C, A_D = areas
```

We inspect the Monte-Carlo draws visually:

```@example propagation
plot(uspec, bnds; size=(500, 600), local_baseline=true)
```

We put in our calculated intensities and their uncertainties (see [`MonteCarloMeasurements.jl` documentation](https://baggepinnen.github.io/MonteCarloMeasurements.jl/stable/api/#MonteCarloMeasurements.:..-Tuple{Any,Any})):

```@example propagation
using MonteCarloMeasurements

I_A = @samples 10_000   1 .. 1.1  # uniform distribution (.. operator)
I_B = @samples 10_000 1.50 ± 0.2  # normal distribution (± operator)
I_C = @samples 10_000 0.75 ± 0.1
I_D = @samples 10_000   1 .. 1.1
nothing # hide
```

!!! warning "Number of particles must match"
    `MonteCarloMeasurements.jl` propagates uncertainty using a Monte-Carlo
    process with a specific number of samples called "particles". By default,
    `MonteCarloMeasurements.jl` produces 2000 particles when defining an
    uncertain number using the operators `±` and `..` (as of version 0.9.5).
    `NoisySignalIntegration` produces 10000 samples by default. You can use the
    [`@samples`](@ref) macro provided by `NoisySignalIntegration.jl` to increase
    the number of samples produced by the `±` and `..` operators as shown in the
    code example above. The number of samples in all uncertain numbers must
    match when calculating combined uncertainties with
    `MonteCarloMeasurements.jl`, so make sure that this condition is met.

Now we can simply calculate with the retrieved areas and defined intensities and
`MonteCarloMeasurements.jl` takes care of the uncertainty calculation:

```@example propagation
R_BA = (A_B/I_B) / (A_A/I_A)
R_CA = (A_C/I_C) / (A_A/I_A)
R_DA = (A_D/I_D) / (A_A/I_A)
nothing # hide
```

This results in the following, displayed as value ± standard deviation:

```@example propagation
R_BA
```

```@example propagation
R_CA
```

```@example propagation
R_DA
```

We can plot the results as histograms to observe the shape of the distributions:

```@example propagation
plot(
    [R_BA.particles R_CA.particles R_DA.particles],
    label=["B:A" "C:A" "D:A"],
    seriestype=:stephist,
    normalize=true,
    xlabel="abundance ratio",
    ylabel="rel. frequency",
    xlim=(0, 10),
    ylim=(0, 2.5),
    fill=true,
    layout=(3, 1)
)
```

Clearly, the resulting distributions are asymmetric and non-Gaussian, so the
standard deviations do not inform about the level of confidence. You can use
[`StatsBase.jl`](https://juliastats.org/StatsBase.jl/stable/) to calculate
percentiles and confidence intervals:

```@example propagation
using StatsBase: percentile

[percentile(R_BA.particles, p) for p in (2.5, 50, 97.5)]
```

A visual representation like the histogram may better convey the range of uncertainty, though.

Alternatives are a box plot:

```@example propagation
using StatsPlots: boxplot

let
    rep(str) = repeat([str], length(R_BA.particles)) 
    x = [rep("B:A") rep("C:A") rep("D:A")]
    y = [R_BA.particles R_CA.particles R_DA.particles]
    boxplot(x, y, ylabel="abundance ratio", label=nothing)
end
```

Or a violin plot:

```@example propagation
using StatsPlots: violin

let
    rep(str) = repeat([str], length(R_BA.particles)) 
    x = [rep("B:A") rep("C:A") rep("D:A")]
    y = [R_BA.particles R_CA.particles R_DA.particles]
    violin(x, y, ylabel="abundance ratio", label=nothing)
end
```
# Package Overview

## Workflow

`NoisySignalIntegration.jl` estimates the uncertainty in numeric integrals based
on the noise level and uncertainty in placing integration bounds. This is
achieved by performing the integration many times while varying the noise and
integration bounds.

The package uses a custom datatype [`Curve`](@ref) to represent the xy-data that
shall be integrated. `Curve` wraps two vectors of identical element type and
length. It was introduced mainly for convenience and simpler plotting.

From a `Curve` object, a [`NoiseSample`](@ref) can be derived. A `NoiseSample`
is required to determine the noise amplitude and autocorrelation length (if the
noise is strongly correlated, as it is often the case in FTIR spectra). With the
noise parameters, a noise model can be constructed. This is either a
[`GaussianNoiseModel`](@ref) (uncorrelated Gaussian noise) or an
[`MvGaussianNoiseModel`](@ref) (correlated Gaussian noise, "Mv" = multivariate).

From the `Curve` object and noise model, an [`UncertainCurve`](@ref) object can
be constructed. It contains random samples of the original data with varying
noise. The `UncertainCurve` object is the first input required for the actual
integration function [`mc_integrate`](@ref).

!!! tip "Crop your data to the relevant region" 
    While your input dataset should contain a somewhat lengthy portion of noise
    for the noise analysis step, you should not include this portion of the data
    in the actual integration, as this will only decrease performance while not
    offering any benefits. You should always [`crop`](@ref) your data to only
    include the relevant signals you want to integrate (see also [Usage
    Guide](@ref crop)).

The second and last input for the integration function is one or several
integration bounds. They are represented by [`UncertainBound`](@ref) objects
which are [defined from one or two probability distributions](@ref bounds_guide)
that encode the uncertainty in the integration bounds.

In the final step, the [`mc_integrate`](@ref) function takes the
`UncertainCurve` and `UncertainBound`(s) and integrates each random sample using
the trapezoidal rule. Each `UncertainBound` yields one uncertain area that is
returned as a
[Particles](https://baggepinnen.github.io/MonteCarloMeasurements.jl/stable/api/#MonteCarloMeasurements.Particles)
object from the
[MonteCarloMeasurements.jl](https://github.com/baggepinnen/MonteCarloMeasurements.jl)
package. [Uncertainty propagation](@ref example_propagation) with the resulting
areas thus works as described in the documentation of
`MonteCarloMeasurements.jl`, merely by performing calculations with `Particles`
objects.

!!! tip "Swapping the core integration function"
    You can swap the core integration function that `mc_integrate` uses to
    something more accurate, if needed. To do so, pass your integration function
    as keyword argument `intfun`. Note, that your integration function needs to
    have the same call signature as [`trapz`](@ref).

    See also: documentation of [`mc_integrate`](@ref).

## Data requirements

Due to the internal workings of the package, the input data needs to fulfill
some basic requirements:

!!! warning "Data must be ordered"
    In general, the data that you analyze must be ordered from low to high
    x-values. If your data is not ordered, you should run the `Base.sort()`
    function on your input data (`Curve` object) once in the beginning of your
    analysis:

```jldoctest
julia> using NoisySignalIntegration

julia> c = Curve([2, 6, 1], [4, 12, 2]);

julia> c = sort(c)
Curve{Int64}, 3 datapoints
(1, 2)
(2, 4)
(6, 12)
```

!!! warning "x-grid must be uniform when analyzing correlated noise"
    Correlated noise can only be analyzed, if the x-data is evenly spaced. If
    this is not the case for your input data, you can use
    [Interpolations.jl](http://juliamath.github.io/Interpolations.jl/latest/) to
    interpolate your data on an even x-grid (see in particular
    [here](http://juliamath.github.io/Interpolations.jl/latest/convenience-construction/)).

Example: Interpolating data on evenly spaced grid

```@example evengrid
using NoisySignalIntegration, Interpolations, Plots

c = NoisySignalIntegration.testdata_3() # test dataset with uneven grid spacing
diff(c.x) # shows the uneven spacing of x-values ↓↓↓
```

```@example evengrid
δx = minimum(diff(c.x))
x_even = collect(minimum(c.x):δx:maximum(c.x))
interp = LinearInterpolation(c.x, c.y)
c_even = Curve(x_even, interp(x_even))

plot(c; label="original data")
plot!(c_even; label="interpolated data")
```