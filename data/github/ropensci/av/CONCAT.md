# av

> R Bindings to FFmpeg

[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![Build Status](https://travis-ci.org/ropensci/av.svg?branch=master)](https://travis-ci.org/ropensci/av)
[![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/github/ropensci/av?branch=master)](https://ci.appveyor.com/project/jeroen/av)
[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/av)](https://cran.r-project.org/package=av)

## Installation

You can install `av` from CRAN

```r
install.packages("av")
```

On Debian/Ubuntu you first need to install [libavfilter-dev](https://packages.debian.org/bullseye/libavfilter-dev)

```
sudo apt-get install libavfilter-dev
```

And on Fedora / CentOS / RHEL you need to install `ffmpeg-devel` from [rpmfusion](https://rpmfusion.org/Configuration). See [instructions here](https://rpmfusion.org/Configuration#Command_Line_Setup_using_rpm) on how to enable rpmfusion via the command line.

```
# Need to enable rpmfusion repository first via link above!
sudo yum install ffmpeg-devel
```

## Demo Video

Generate a demo video with some random plots and free [demo music](https://freemusicarchive.org/music/Synapsis/~/Wonderland):

```r
av::av_demo()
```

This demo is totally lame, please open a PR with something better (in base R!).

## Using gganimate

You can use `av_encode_video()` as the renderer in gganimate:

```r
# Create the gganimate plot
library(gganimate)
library(transformr)
p <- ggplot(airquality, aes(Day, Temp)) + 
  geom_line(size = 2, colour = 'steelblue') + 
  transition_states(Month, 4, 1) + 
  shadow_mark(size = 1, colour = 'grey')

# Render and show the video
q <- 2
df <- animate(p, renderer = av_renderer('animation.mp4'), width = 720*q, height = 480*q, res = 72*q, fps = 25)
utils::browseURL('animation.mp4')
```

## Video Filters

You can add a custom [ffmpeg video filter chain](https://ffmpeg.org/ffmpeg-filters.html#Video-Filters). For example this will negate the colors, and applies an orange fade-in effect to the first 15 frames.

```r
# Continue on the example above
filter_render <- av_renderer('orange.mp4', vfilter = 'negate=1, fade=in:0:15:color=orange')
df <- animate(p, renderer = filter_render, width = 720*q, height = 480*q, res = 72*q, fps = 25)
av::av_media_info('orange.mp4')
utils::browseURL('orange.mp4')
```

Filters can also affect the final fps of the video. For example this filter will double fps because it halves presentation the timestamp (pts) of each frame. Hence the output framerate is actually 50!

```r
fast_render <- av_renderer('fast.mp4', vfilter = "setpts=0.5*PTS")
df <- animate(p, renderer = fast_render, fps = 25)
av::av_media_info('fast.mp4')
utils::browseURL('fast.mp4')
```

## Capture Graphics (without gganimate)

Instead of using gganimate, we can use `av_capture_graphics()` to automatically record R graphics and turn them into a video. This example makes 12 plots and adds an interpolation filter to smoothen the transitions between the frames.

```r
library(gapminder)
library(ggplot2)
makeplot <- function(){
  datalist <- split(gapminder, gapminder$year)
  lapply(datalist, function(data){
    p <- ggplot(data, aes(gdpPercap, lifeExp, size = pop, color = continent)) +
      scale_size("population", limits = range(gapminder$pop)) + geom_point() + ylim(20, 90) +
      scale_x_log10(limits = range(gapminder$gdpPercap)) + ggtitle(data$year) + theme_classic()
    print(p)
  })
}

# Play 1 plot per sec, and use an interpolation filter to convert into 10 fps
video_file <- file.path(tempdir(), 'output.mp4')
av::av_capture_graphics(makeplot(), video_file, 1280, 720, res = 144, vfilter = 'framerate=fps=10')
av::av_media_info(video_file)
utils::browseURL(video_file)
```
---
title: "Spectrograms in R using the 'av' package"
author: "Jeroen Ooms"
date: "`r Sys.Date()`"
---

```{r setup, include = FALSE}
library(av)
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
```

Calculate the frequency data and plot the spectrogram:

```{r plot_spectrogram}
# Demo sound included with av
wonderland <- system.file('samples/Synapsis-Wonderland.mp3', package='av')

# Read first 5 sec of demo
fft_data <- read_audio_fft(wonderland, end_time = 5.0)
plot(fft_data)
```

You can turn off dark mode to use the default R colors:

```{r plot_spectrogram_white}
plot(fft_data, dark = FALSE)
```

## Spectrogram video

You can also create a spectrogram video like this:

```{r spectrogram_video}
# Create new audio file with first 5 sec
av_audio_convert(wonderland, 'short.mp3', total_time = 5)
av_spectrogram_video('short.mp3', output = 'spectrogram.mp4', width = 1280, height = 720, res = 144)
```

<video width="100%" controls>
<source src="spectrogram.mp4" type="video/mp4">
Your browser does not support the video tag.
</video>

## Compare with tuneR/signal

For comparison, we show how the same thing can be achieved with the `tuneR` package:

```{r convert_to_wav}
# Read wav with tuneR
data <- tuneR::readMP3('short.mp3')

# demean to remove DC offset
snd <- data@left - mean(data@left)
```

We then use the signal package to calculate the spectrogram with similar parameters as av:

```{r tuner_specgram}
# create spectrogram
spec <- signal::specgram(x = snd, n = 1024, Fs = data@samp.rate, overlap = 1024 * 0.75)

# normalize and rescale to dB
P <- abs(spec$S)
P <- P/max(P)

out <- pmax(1e-6, P)
dim(out) <- dim(P)
out <- log10(out) / log10(1e-6)

# plot spectrogram
image(x = spec$t, y = spec$f, z = t(out), ylab = 'Freq [Hz]', xlab = 'Time [s]', useRaster=TRUE)
```

## Compare with seewave

Compare spectrograms using the `tico` audio sample included with the seewave package:

```{r seewave_spectrogram}
library(seewave)
library(ggplot2)
data(tico)
ggspectro(tico, ovlp = 50) + geom_tile(aes(fill = amplitude)) + stat_contour()
```

To use av, we first save the wav file and then create spectrogram:

```{r av_tico}
# First export wav file
savewav(tico, filename = 'tico.wav')
plot(read_audio_fft('tico.wav'))
```

## Compare with phonTools

Use the audio sample included with phonTools:

```{r phontools_spectrogram}
library(phonTools)
data(sound)
spectrogram(sound, maxfreq = sound$fs/2)
```

Save the wav file and then create spectrogram. We match the default window function from phonTools:

```{r av_phontools}
phonTools::writesound(sound, 'sound.wav')
plot(read_audio_fft('sound.wav', window = phonTools::windowfunc(1024, 'kaiser')))
```
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/capture.R
\name{capturing}
\alias{capturing}
\alias{av_capture_graphics}
\alias{av_spectrogram_video}
\title{Record Video from Graphics Device}
\usage{
av_capture_graphics(
  expr,
  output = "output.mp4",
  width = 720,
  height = 480,
  framerate = 1,
  vfilter = "null",
  audio = NULL,
  verbose = TRUE,
  ...
)

av_spectrogram_video(
  audio,
  output = "output.mp4",
  framerate = 25,
  verbose = TRUE,
  ...
)
}
\arguments{
\item{expr}{an R expression that generates the graphics to capture}

\item{output}{name of the output file. File extension must correspond to a known
container format such as \code{mp4}, \code{mkv}, \code{mov}, or \code{flv}.}

\item{width}{width in pixels of the graphics device}

\item{height}{height in pixels of the graphics device}

\item{framerate}{video framerate in frames per seconds. This is the input fps, the
output fps may be different if you specify a filter that modifies speed or interpolates
frames.}

\item{vfilter}{a string defining an ffmpeg filter graph. This is the same parameter
as the \code{-vf} argument in the \code{ffmpeg} command line utility.}

\item{audio}{path to media file with audio stream}

\item{verbose}{emit some output and a progress meter counting processed images. Must
be \code{TRUE} or \code{FALSE} or an integer with a valid \link{av_log_level}.}

\item{...}{extra graphics parameters passed to \code{\link[=png]{png()}}}
}
\description{
Runs the expression and captures all plots into a video. The \link{av_spectrogram_video}
function is a wrapper that plots data from \link{read_audio_fft} with a moving bar and
background audio.
}
\examples{
\donttest{
library(gapminder)
library(ggplot2)
makeplot <- function(){
  datalist <- split(gapminder, gapminder$year)
  lapply(datalist, function(data){
    p <- ggplot(data, aes(gdpPercap, lifeExp, size = pop, color = continent)) +
      scale_size("population", limits = range(gapminder$pop)) + geom_point() + ylim(20, 90) +
      scale_x_log10(limits = range(gapminder$gdpPercap)) + ggtitle(data$year) + theme_classic()
    print(p)
  })
}

# Play 1 plot per sec, and use an interpolation filter to convert into 10 fps
video_file <- file.path(tempdir(), 'output.mp4')
av_capture_graphics(makeplot(), video_file, 1280, 720, res = 144, vfilter = 'framerate=fps=10')
av::av_media_info(video_file)
# utils::browseURL(video_file)}
}
\seealso{
Other av: 
\code{\link{demo}()},
\code{\link{encoding}},
\code{\link{formats}},
\code{\link{info}},
\code{\link{logging}},
\code{\link{read_audio_fft}()}
}
\concept{av}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fft.R
\name{read_audio_fft}
\alias{read_audio_fft}
\alias{read_audio_bin}
\title{Read audio binary and frequency data}
\usage{
read_audio_fft(
  audio,
  window = hanning(1024),
  overlap = 0.75,
  sample_rate = NULL,
  start_time = NULL,
  end_time = NULL
)

read_audio_bin(
  audio,
  channels = NULL,
  sample_rate = NULL,
  start_time = NULL,
  end_time = NULL
)
}
\arguments{
\item{audio}{path to the input sound or video file containing the audio stream}

\item{window}{vector with weights defining the moving \link[=hanning]{fft window function}.
The length of this vector is the size of the window and hence determines the output
frequency range.}

\item{overlap}{value between 0 and 1 of overlap proportion between moving fft windows}

\item{sample_rate}{downsample audio to reduce FFT output size. Default keeps sample
rate from the input file.}

\item{start_time, end_time}{position (in seconds) to cut input stream to be processed.}

\item{channels}{number of output channels, set to 1 to convert to mono sound}
}
\description{
Reads raw audio data from any common audio or video format. Use \link{read_audio_bin} to
get raw PCM audio samples, or \link{read_audio_fft} to stream-convert directly into
frequency domain (spectrum) data using FFmpeg built-in FFT.
}
\details{
Currently \link{read_audio_fft} automatically converts input audio to mono channel such
that we get a single matrix. Use the \code{plot()} method on data returned by \link{read_audio_fft}
to show the spectrogram. The \link{av_spectrogram_video} generates a video that plays
the audio while showing an animated spectrogram with moving status bar, which is
very cool.
}
\examples{
# Use a 5 sec fragment
wonderland <- system.file('samples/Synapsis-Wonderland.mp3', package='av')

# Read initial 5 sec as as frequency spectrum
fft_data <- read_audio_fft(wonderland, end_time = 5.0)
dim(fft_data)

# Plot the spectrogram
plot(fft_data)

# Show other parameters
dim(read_audio_fft(wonderland, end_time = 5.0, hamming(2048)))
dim(read_audio_fft(wonderland, end_time = 5.0, hamming(4096)))
}
\seealso{
Other av: 
\code{\link{capturing}},
\code{\link{demo}()},
\code{\link{encoding}},
\code{\link{formats}},
\code{\link{info}},
\code{\link{logging}}
}
\concept{av}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/formats.R
\name{formats}
\alias{formats}
\alias{av_encoders}
\alias{av_decoders}
\alias{av_filters}
\alias{av_muxers}
\alias{av_demuxers}
\title{AV Formats}
\usage{
av_encoders()

av_decoders()

av_filters()

av_muxers()

av_demuxers()
}
\description{
List supported filters, codecs and container formats.
}
\details{
Encoders and decoders convert between raw video/audio frames and compressed stream
data for storage or transfer. However such a compressed data stream by itself does
not constitute a valid video format yet. Muxers are needed to interleave one or more
audio/video/subtitle streams, along with timestamps, metadata, etc, into a proper
file format, such as mp4 or mkv.

Conversely, demuxers are needed to read a file format into the seperate data streams
for subsequent decoding into raw audio/video frames. Most operating systems natively
support demuxing and decoding common formats and codecs, needed to play those videos.
However for encoding and muxing such videos, ffmpeg must have been configured with
specific external libraries for a given codec or format.
}
\seealso{
Other av: 
\code{\link{capturing}},
\code{\link{demo}()},
\code{\link{encoding}},
\code{\link{info}},
\code{\link{logging}},
\code{\link{read_audio_fft}()}
}
\concept{av}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/images.R
\name{av_video_images}
\alias{av_video_images}
\title{Convert video to images}
\usage{
av_video_images(video, destdir = tempfile(), format = "jpg", fps = NULL)
}
\arguments{
\item{video}{an input video}

\item{destdir}{directory where to save the png files}

\item{format}{image format such as \code{png} or \code{jpeg}, must be available from \code{av_encoders()}}

\item{fps}{sample rate of images. Use \code{NULL} to get all images.}
}
\description{
Splits a video file in a set of image files. Default image format is
jpeg which has good speed and compression. Use \code{format = "png"} for
losless images.
}
\details{
For large input videos you can set fps to sample only a limited number
of images per second. This also works with fractions, for example \code{fps = 0.2}
will output one image for every 5 sec of video.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/encode.R
\name{encoding}
\alias{encoding}
\alias{av_encode_video}
\alias{av}
\alias{av_video_convert}
\alias{av_audio_convert}
\title{Encode or Convert Audio / Video}
\usage{
av_encode_video(
  input,
  output = "output.mp4",
  framerate = 24,
  vfilter = "null",
  codec = NULL,
  audio = NULL,
  verbose = TRUE
)

av_video_convert(video, output = "output.mp4", verbose = TRUE)

av_audio_convert(
  audio,
  output = "output.mp3",
  format = NULL,
  channels = NULL,
  sample_rate = NULL,
  start_time = NULL,
  total_time = NULL,
  verbose = TRUE
)
}
\arguments{
\item{input}{a vector with image or video files. A video input file is treated
as a series of images. All input files should have the same width and height.}

\item{output}{name of the output file. File extension must correspond to a known
container format such as \code{mp4}, \code{mkv}, \code{mov}, or \code{flv}.}

\item{framerate}{video framerate in frames per seconds. This is the input fps, the
output fps may be different if you specify a filter that modifies speed or interpolates
frames.}

\item{vfilter}{a string defining an ffmpeg filter graph. This is the same parameter
as the \code{-vf} argument in the \code{ffmpeg} command line utility.}

\item{codec}{name of the video codec as listed in \link{av_encoders}. The
default is \code{libx264} for most formats, which usually the best choice.}

\item{audio}{audio or video input file with sound for the output video}

\item{verbose}{emit some output and a progress meter counting processed images. Must
be \code{TRUE} or \code{FALSE} or an integer with a valid \link{av_log_level}.}

\item{video}{input video file with optionally also an audio track}

\item{format}{a valid format name from the list of \code{av_muxers()}. Default
\code{NULL} tries to guess a format from the file extension.}

\item{channels}{number of output channels. Default \code{NULL} is to match input}

\item{sample_rate}{output sampling rate. Default \code{NULL} is to match input}

\item{start_time}{number greater than 0, seeks in the input file to position.}

\item{total_time}{approximate number of seconds at which to limit the duration
of the output file.}
}
\description{
Encodes a set of images into a video, using custom container format, codec, fps,
\href{https://ffmpeg.org/ffmpeg-filters.html#Video-Filters}{video filters}, and audio
track. If input contains video files, this effectively combines and converts them
to the specified output format.
}
\details{
The target container format is automatically determined from the file extension of
the output file, for example \code{mp4}, \code{mkv}, \code{mov}, or \code{flv}. Most systems also
support \code{gif} output, but the compression~quality for gif is quite bad.
The \href{https://cran.r-project.org/package=gifski}{gifski} package is better suited
for generating animated gif files.

It is recommended to use let ffmpeg choose the suitable codec for a given container
format. Most video formats default to the \code{libx264} video codec which has excellent
compression and works on all modern \href{https://caniuse.com/#search=h264}{browsers},
operating systems, and digital TVs.

It is safe to interrupt the encoding process by pressing CTRL+C, or via \link{setTimeLimit}.
When the encoding is interrupted, the output stream is properly finalized and all open
files and resources are properly closed.
}
\seealso{
Other av: 
\code{\link{capturing}},
\code{\link{demo}()},
\code{\link{formats}},
\code{\link{info}},
\code{\link{logging}},
\code{\link{read_audio_fft}()}
}
\concept{av}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/demo.R
\name{demo}
\alias{demo}
\alias{av_demo}
\title{Demo Video}
\usage{
av_demo(
  output = "demo.mp4",
  width = 960,
  height = 720,
  framerate = 5,
  verbose = TRUE,
  ...
)
}
\arguments{
\item{output}{name of the output file. File extension must correspond to a known
container format such as \code{mp4}, \code{mkv}, \code{mov}, or \code{flv}.}

\item{width}{width in pixels of the graphics device}

\item{height}{height in pixels of the graphics device}

\item{framerate}{video framerate in frames per seconds. This is the input fps, the
output fps may be different if you specify a filter that modifies speed or interpolates
frames.}

\item{verbose}{emit some output and a progress meter counting processed images. Must
be \code{TRUE} or \code{FALSE} or an integer with a valid \link{av_log_level}.}

\item{...}{other parameters passed to \link{av_capture_graphics}.}
}
\description{
Generates random video for testing purposes.
}
\seealso{
Other av: 
\code{\link{capturing}},
\code{\link{encoding}},
\code{\link{formats}},
\code{\link{info}},
\code{\link{logging}},
\code{\link{read_audio_fft}()}
}
\concept{av}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/logging.R
\name{logging}
\alias{logging}
\alias{av_log_level}
\title{Logging}
\usage{
av_log_level(set = NULL)
}
\arguments{
\item{set}{new \href{https://www.ffmpeg.org/doxygen/4.0/group__lavu__log__constants.html}{log level} value}
}
\description{
Get or set the \href{https://www.ffmpeg.org/doxygen/4.0/group__lavu__log__constants.html}{log level}.
}
\seealso{
Other av: 
\code{\link{capturing}},
\code{\link{demo}()},
\code{\link{encoding}},
\code{\link{formats}},
\code{\link{info}},
\code{\link{read_audio_fft}()}
}
\concept{av}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/winfunc.R
\name{window functions}
\alias{window functions}
\alias{hanning}
\alias{hamming}
\alias{blackman}
\alias{bartlett}
\alias{welch}
\alias{flattop}
\alias{bharris}
\alias{bnuttall}
\alias{sine}
\alias{nuttall}
\alias{bhann}
\alias{lanczos}
\alias{gauss}
\alias{tukey}
\alias{dolph}
\alias{cauchy}
\alias{parzen}
\alias{bohman}
\title{Window functions}
\usage{
hanning(n)

hamming(n)

blackman(n)

bartlett(n)

welch(n)

flattop(n)

bharris(n)

bnuttall(n)

sine(n)

nuttall(n)

bhann(n)

lanczos(n)

gauss(n)

tukey(n)

dolph(n)

cauchy(n)

parzen(n)

bohman(n)
}
\arguments{
\item{n}{size of the window (number of weights to generate)}
}
\description{
Several common \href{https://en.wikipedia.org/wiki/Window_function}{windows function}
generators. The functions return a vector of weights to use in \link{read_audio_fft}.
}
\examples{
# Window functions
plot(hanning(1024), type = 'l', xlab = 'window', ylab = 'weight')
lines(hamming(1024), type = 'l', col = 'red')
lines(bartlett(1024), type = 'l', col = 'blue')
lines(welch(1024), type = 'l', col = 'purple')
lines(flattop(1024), type = 'l', col = 'darkgreen')
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/info.R
\name{info}
\alias{info}
\alias{av_media_info}
\alias{av_video_info}
\title{Video Info}
\usage{
av_media_info(file)
}
\arguments{
\item{file}{path to an existing file}
}
\description{
Get video info such as width, height, format, duration and framerate.
This may also be used for audio input files.
}
\seealso{
Other av: 
\code{\link{capturing}},
\code{\link{demo}()},
\code{\link{encoding}},
\code{\link{formats}},
\code{\link{logging}},
\code{\link{read_audio_fft}()}
}
\concept{av}
