# Package pamlr

This package combines a set of functions for analysing multisensor geolocator data such as Pressure, Activity, Magnetism, Temperature or Light. This includes functions for importing and visualising any or all of these sensor data, formatting the data for analysis (with some functions already setup for birds), deriving classifications (using cluster analysis or hidden markove models), and finally for comparing classification accuracy between different models.

## Access the package Manual 

The manual is the best way to get started with the package, and can be found at the following link:
[https://kiranlda.github.io/PAMLrManual/index.html](https://kiranlda.github.io/PAMLrManual/index.html)


## Installation

### Prerequisites

To install this package from github, make sure the user first have `devtools` installed.

```r
install.packages("devtools")
```
Then the package can be installed from github:

```r
devtools::install_github("KiranLDA/pamlr")
```
## Import data

Data can be imported as follows:

```r
PAM_data = create_import(pathname = "C:/Put/your/path/here",
                     measurements = c(".pressure", 
                                      ".glf",
                                      ".acceleration", 
                                      ".temperature", 
                                      ".magnetic")`
```

An example dataset is 

```r
data(hoopoe)
```

## Cropping the data

Note that very often, a logger continues to record data before and after it is removed from a bird. For example, it might be transported in a rucksack or left in a laboratory until data are downloaded. It's therefore important to remove these incorrect datapoints. This can be done using `create_crop`.

```r
# make sure the cropping period is in the correct date format
start = as.POSIXct("2016-07-01","%Y-%m-%d", tz="UTC")
end = as.POSIXct("2017-06-01","%Y-%m-%d", tz="UTC")

# Crop the data
PAM_data= create_crop(hoopoe,start,end)
```

## Visualising data

### Quick multiplots

For a quick look at the data, it's possible to use `plot_timeseries`. The user can specify which arguments to use, using `measurements`. There's a choice between different combinations of `"pressure"`, `"light"`, `"acceleration"`, `"temperature"` and `"magnetic"`. You can add any parameters from `?plot`, here I illustrate it with `col="cornflowerblue"` and by showing how to  restrict the x-axis limits `xlim` with the date format, to zoom into the post breeding mihration period of a hoopoe

```r
par(mar=c(3,4,0.5,0.5))
plot_timeseries(hoopoe, col="cornflowerblue",
                measurements = c("pressure", "light", "acceleration"),
                xlim=c(as.POSIXct("2016-08-20","%Y-%m-%d", tz="UTC"),
                        as.POSIXct("2016-09-01","%Y-%m-%d", tz="UTC")))
```
<img align="center" src="https://github.com/KiranLDA/PAMLrManual/blob/master/_bookdown_files/PAMLrManual_files/figure-html/unnamed-chunk-14-1.png">


### Interactive timeseries

To have a better overview of the data, it is possible to create interactive `plot_interactive_timeseries()` plots which allow the user to compare different measurements recorded by the logger. These might for instance include light, temperature, pressure, activity, pitch and magnetism. 

If you are **working from Rstudio**, this bit of code should be run:

```r
# In Rstudio, it will display in the viewer by default and use a lot of ram, and is better in html
backup_options <- options() 
options(viewer=NULL) # ensure it is viewed in internet browser
plot_interactive_timeseries(dta = PAM_data) # plot
options(backup_options) # restore previous viewer settings
```

If you are **working from base R** use this instead:

To save space here we only plot only one variable - pressure .

```r
plot_interactive_timeseries(dta = PAM_data, toPLOT = c("pressure")) 
```

The reason there is additional code for Rstudio, is that by default it will open this graphic in the viewer pane and use up a lot of ram. This  additional code allows the user to open this window in a browser instead of r studio, and the file can later be saved as an html file.

With this interactive plot, the user can then zoom in and out of different plots to help get a feel for the data. For instance, this is a great way of seeing changes in the data which might be due to a logger being in a rucksack and no longer on the birds, or to look at how acticity or pressure might look during migration periods.

It is possible to select areas to zoom into by right clicking and highighting certain regions, and to double click to zoom out. All plots are synched to the same time period and have a timeline at the bottom to increase or decrease the time over which the data is observed.

### Sensor images

Actograms are often used to plot activity over time at different hours of the day. However, the same approach can be used to plot any sensor data, not just activity. For simplicity, we name these “sensor images”. Plotting all sensors side by side is an important step for visualising data and developing an understanding of data patterns, and to start thinking about the behaviours that may be driving the observed patterns. __PAMLr__ offers a function `plot_sensorimage()`for plotting sensor images, which can be implemented as follows.

```r
# Create plots with 3 together (mfrow)
par( mfrow= c(1,3), oma=c(0,2,0,6))

par(mar =  c(4,2,4,2))
plot_sensorimage(PAM_data$acceleration$date, ploty=FALSE,
          PAM_data$acceleration$act, main = "Activity",
          col=c("black",viridis::cividis(90)), cex=1.2, cex.main = 2)

par(mar =  c(4,2,4,2))
plot_sensorimage(PAM_data$pressure$date, plotx=TRUE, ploty=FALSE, labely=FALSE,
          PAM_data$pressure$obs,  main="Pressure",
          col=c("black",viridis::cividis(90)), cex=1.2, cex.main = 2)

par(mar =  c(4,2,4,2))
plot_sensorimage(PAM_data$temperature$date, labely=FALSE,
          PAM_data$temperature$obs,  main="Temperature",
          col=c("black",viridis::cividis(90)), cex=1.2, cex.main = 2)


```
<img align="center" src="https://github.com/KiranLDA/PAMLrManual/blob/master/_bookdown_files/PAMLrManual_files/figure-html/unnamed-chunk-17-1.png">

### Histograms and three-dimensional scatterplots

Histograms can provide a first impression of whether some of the data may be aggregated and therefore clustered. Indeed, sensor images may not always well-suited for visualising tri-axial data, such as magnetic field or acceleration. By plotting data in three dimensions (hereafter “3D”) using the function `plot_interactive_3d` it’s possible to find patterns or clusters of datapoints which would not otherwise be apparent in the data. Here we plot magnetic data. 

```r
plot_interactive_3d(PAM_data$magnetic$gX, PAM_data$magnetic$gY, PAM_data$magnetic$gZ,
       xlab= "X", ylab= "Y", zlab= "Z",
       xlim=c(-3000,3000), ylim=c(-3000,3000), zlim=c(-3000,3000))
```

### Spherical projections

#### g-sphere

A _g-sphere_ is a method of visualising tri-axial *acceleration* data. This involves centering the data and plotting it on a sphere.The function `calculate_triaxial_acceleration` allows the user to center this data (as well as calculating  pitch, roll and yaw from the data)

```r
# plot an g-phere
calibration = calculate_triaxial_acceleration(dta = PAM_data$magnetic)
plot_interactive_sphere(x = calibration$centered_accx,
          y = calibration$centered_accy,
          z = calibration$centered_accz,
          ptcol = "royalblue",
          ptsize = 0.03,
          linecolor ="orange",
          spherecolor="orange",
          arrows=TRUE)


```

#### m-sphere

An _m-sphere_ is a method of visualising tri-axial *magnetometer* data. This involves centering the data and correcting the data, before  plotting it on a sphere.The function `calculate_triaxial_magnetic` calibrates the data. This provides the animal's bearing.

```r
# plot a m-phere
calibration = calculate_triaxial_magnetic(dta = PAM_data$magnetic)
plot_interactive_sphere(x = calibration$calib_magx,
          y = calibration$calib_magy,
          z = calibration$calib_magz,
          ptcol = "orange",
          ptsize = 0.03,
          linecolor ="royalblue",
          spherecolor="royalblue",
          arrows=TRUE,
          cex=2)
```


## Classifying migratory flapping flight in Hoopoes  (an example classification)

## Performing the classification

Because flapping is widespread in birds, **PAMLr** integrates a pre-defined function `classify_flap()` to classify this behaviour.

This function assumes that if the bird has displayed **high activity** for x number of consecutive minutes, then it is flapping. It is therefore important to think about what constitutes high activity and how long this period should be. At the moment, the function uses **k-means clustering** to identify the threshold between high and low activity. Using `toPLOT = TRUE` then allows the user to see where that threshold was drawn. The period of high activity is set by default to `period = 3`. This is because activity is recorded (on this logger) every 5 minutes and we assume that after an hour of high activity, the bird must be flapping. 

Thus "high activity duration" / "data resolution" = "period" and 60 minutes / 5 minutes = period of 12.

```r
# Classify behaviour
behaviour = classify_flap(dta = PAM_data$acceleration, period = 12)
str(behaviour)
```
<img align="center" src="https://github.com/KiranLDA/PAMLrManual/blob/master/_bookdown_files/PAMLrManual_files/figure-html/unnamed-chunk-50-1.png">


This classification therefore provides different pieces of infomration.

* **timetable** shows when a migratory flapping flight started and stopped, and how long it lasted (in hours)
* **classification** is the output from the classification. In this case, each cluster/classs/state is represented by numbers between one 1 and 4. To find out what behaviour each of these numbers represent, we can refer to **low_movement**, **high_movement**, **migration** and **no_movement** 
* **threshold** represents the threshold between high and low activity.

Using these information, it's therefore possible to plot the classification:


```r
#Plot behaviour
col=col=c("black","royalblue4","brown","gold")
index= 7300:8000
plot(PAM_data$acceleration$date[index],PAM_data$acceleration$act[index],
  type="l", xlab="Date", ylab="Activity")
points(PAM_data$acceleration$date[index],PAM_data$acceleration$act[index],
  col=col[behaviour$classification+1][index], 
  pch=16,)
legend( PAM_data$acceleration$date[index[1]],60 , 
        c("No activity", "Low activity", "High activity", "Migration" ) ,
        col = col[c(behaviour$no_movement, behaviour$low_movement,
                    behaviour$high_movement, behaviour$migration)+1],
        pch=20)
```
<img align="center" src="https://github.com/KiranLDA/PAMLrManual/blob/master/docs/08-flapping_files/figure-html/unnamed-chunk-5-1.png">

## Plot the classification as a sensor image

Another way of looking at a classification is to use a sensor image of the results and to plot it side by side with the raw data to see if the same patterns are being picked out. We can also add (for instance sunset and sunrise events)

```r
par(mfrow= c(1,3), # number of panels
    oma=c(0,2,0,6), # outer margin around all panels
    mar =  c(4,1,4,1)) # inner margin around individual fivure

plot_sensorimage(PAM_data$acceleration$date, ploty=FALSE,
          PAM_data$acceleration$act, main = "Activity",
          col=c("black",viridis::cividis(90)), cex=1.2, cex.main = 2)
legend("bottomright",cex=1.2,
   c("No Activity", "Low Activity", "High Activity" ) , fill = c("black","royalblue3", "orange"), xpd = NA)


plot_sensorimage(PAM_data$pressure$date, plotx=TRUE, ploty=FALSE, labely=FALSE,
          PAM_data$pressure$obs,  main="Pressure",
          col=c("black",viridis::cividis(90)), cex=1.2, cex.main = 2)
legend("bottomright",cex=1.2,
   c("Low Pressure", "High Pressure" ) , fill = c("royalblue3", "orange"), xpd = NA)


plot_sensorimage(PAM_data$acceleration$date, labely=FALSE,
          behaviour$classification, 
          main="Classification",
          col=col,
          cex=1.2, cex.main = 2)
legend("bottomright",cex=1.2,
  # grconvertX(1, "device"), grconvertY(1, "device"),
   c("Resting", "Active", "Flapping", "Migrating" ) , fill = col, xpd = NA)


```

<img align="center" src="https://github.com/KiranLDA/PAMLrManual/blob/master/docs/08-flapping_files/figure-html/unnamed-chunk-6-1.png">



## Authors

Kiran Dhanjal-Adams

## Relevant references

Dhanjal-Adams, K.L., Bauer, S., Emmenegger, T., Hahn, S., Lisovski, S., & Liechti, F. (2018) [Spatiotemporal Group Dynamics in a Long-Distance Migratory Bird](https://www.cell.com/current-biology/fulltext/S0960-9822(18)30845-5?_returnURL=https%3A%2F%2Flinkinghub.elsevier.com%2Fretrieve%2Fpii%2FS0960982218308455%3Fshowall%3Dtrue). Current Biology, 28, 2824–2830.e3. 

Liechti, F., Bauer, S., Dhanjal-Adams, K.L., Emmenegger, T., Zehtindjiev, P., & Hahn, S. (2018) [Miniaturized multi-sensor loggers provide new insight into year-round flight behaviour of small trans-Sahara avian migrants](https://movementecologyjournal.biomedcentral.com/articles/10.1186/s40462-018-0137-1). Movement Ecology, 6, 19. 

## License

This project is licensed under the GNU General Public License version 3 - see the [LICENSE](https://github.com/KiranLDA/PAMLr/blob/master/LICENSE) file for details
---
title: "pamlr vignette"
author: "Kiran Dhanjal-Adams"
date: "`r Sys.Date()`"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{pamlr vignette}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---



# Package pamlr

This package combines a set of functions for analysing multisensor geolocator data such as Pressure, Activity, Magnetism, Temperature or Light. This includes functions for importing and visualising any or all of these sensor data, formatting the data for analysis (with some functions already setup for birds), deriving classifications (using cluster analysis or hidden markove models), and finally for comparing classification accuracy between different models.

## Access the package Manual 

The manual is the best way to get started with the package, and can be found at the following link:
[https://kiranlda.github.io/PAMLrManual/index.html](https://kiranlda.github.io/PAMLrManual/index.html)


## Installation

### Prerequisites

To install this package from github, make sure the user first have `devtools` installed.

```r
install.packages("devtools")
```
Then the package can be installed from github:

```r
devtools::install_github("KiranLDA/PAMLr")
```
## Import data

Data can be imported as follows:

```r
library(PAMLr)
PAM_data = importPAM(pathname = "C:/Put/your/path/here",
                     measurements = c(".pressure", 
                                      ".glf",
                                      ".acceleration", 
                                      ".temperature", 
                                      ".magnetic")`
```

An example dataset is 

```r
data(hoopoe)
```

## Cropping the data

Note that very often, a logger continues to record data before and after it is removed from a bird. For example, it might be transported in a rucksack or left in a laboratory until data are downloaded. It's therefore important to remove these incorrect datapoints. This can be done using `cutPAM`.

```r
# make sure the cropping period is in the correct date format
start = as.POSIXct("2016-07-01","%Y-%m-%d", tz="UTC")
end = as.POSIXct("2017-06-01","%Y-%m-%d", tz="UTC")

# Crop the data
PAM_data= cutPAM(hoopoe,start,end)
```

## Visualising data

### Quick multiplots

For a quick look at the data, it's possible to use `quickPLOT`. The user can specify which arguments to use, using `measurements`. There's a choice between different combinations of `"pressure"`, `"light"`, `"acceleration"`, `"temperature"` and `"magnetic"`. You can add any parameters from `?plot`, here I illustrate it with `col="cornflowerblue"` and by showing how to  restrict the x-axis limits `xlim` with the date format, to zoom into the post breeding mihration period of a hoopoe

```r
par(mar=c(3,4,0.5,0.5))
quickPLOT(hoopoe, col="cornflowerblue",
          measurements = c("pressure", "light", "acceleration"),
          xlim=c(as.POSIXct("2016-08-20","%Y-%m-%d", tz="UTC"),
                 as.POSIXct("2016-09-01","%Y-%m-%d", tz="UTC")))
```
<img align="center" src="https://github.com/KiranLDA/PAMLrManual/blob/master/_bookdown_files/PAMLrManual_files/figure-html/unnamed-chunk-14-1.png">


### Interactive timeseries

To have a better overview of the data, it is possible to create interactive `dygraphPAM()` plots which allow the user to compare different measurements recorded by the logger. These might for instance include light, temperature, pressure, activity, pitch and magnetism. 

If you are **working from Rstudio**, this bit of code should be run:

```r
# In Rstudio, it will display in the viewer by default and use a lot of ram, and is better in html
backup_options <- options() 
options(viewer=NULL) # ensure it is viewed in internet browser
dygraphPAM(dta = PAM_data) # plot
options(backup_options) # restore previous viewer settings
```

If you are **working from base R** use this instead:

To save space here we only plot only one variable - pressure .

```r
dygraphPAM(dta = PAM_data, toPLOT = c("pressure")) 
```

The reason there is additional code for Rstudio, is that by default it will open this graphic in the viewer pane and use up a lot of ram. This  additional code allows the user to open this window in a browser instead of r studio, and the file can later be saved as an html file.

With this interactive plot, the user can then zoom in and out of different plots to help get a feel for the data. For instance, this is a great way of seeing changes in the data which might be due to a logger being in a rucksack and no longer on the birds, or to look at how acticity or pressure might look during migration periods.

It is possible to select areas to zoom into by right clicking and highighting certain regions, and to double click to zoom out. All plots are synched to the same time period and have a timeline at the bottom to increase or decrease the time over which the data is observed.

### Sensor images

Actograms are often used to plot activity over time at different hours of the day. However, the same approach can be used to plot any sensor data, not just activity. For simplicity, we name these “sensor images”. Plotting all sensors side by side is an important step for visualising data and developing an understanding of data patterns, and to start thinking about the behaviours that may be driving the observed patterns. __PAMLr__ offers a function `sensorIMG()`for plotting sensor images, which can be implemented as follows.

```r
# Create plots with 3 together (mfrow)
par( mfrow= c(1,3), oma=c(0,2,0,6))

par(mar =  c(4,2,4,2))
sensorIMG(PAM_data$acceleration$date, ploty=FALSE,
          PAM_data$acceleration$act, main = "Activity",
          col=c("black",viridis::cividis(90)), cex=1.2, cex.main = 2)

par(mar =  c(4,2,4,2))
sensorIMG(PAM_data$pressure$date, plotx=TRUE, ploty=FALSE, labely=FALSE,
          PAM_data$pressure$obs,  main="Pressure",
          col=c("black",viridis::cividis(90)), cex=1.2, cex.main = 2)

par(mar =  c(4,2,4,2))
sensorIMG(PAM_data$temperature$date, labely=FALSE,
          PAM_data$temperature$obs,  main="Temperature",
          col=c("black",viridis::cividis(90)), cex=1.2, cex.main = 2)


```
<img align="center" src="https://github.com/KiranLDA/PAMLrManual/blob/master/_bookdown_files/PAMLrManual_files/figure-html/unnamed-chunk-17-1.png">

### Histograms and three-dimensional scatterplots

Histograms can provide a first impression of whether some of the data may be aggregated and therefore clustered. Indeed, sensor images may not always well-suited for visualising tri-axial data, such as magnetic field or acceleration. By plotting data in three dimensions (hereafter “3D”) using the function `pam3D` it’s possible to find patterns or clusters of datapoints which would not otherwise be apparent in the data. Here we plot magnetic data. 

```r
pam3D(PAM_data$magnetic$gX, PAM_data$magnetic$gY, PAM_data$magnetic$gZ,
       xlab= "X", ylab= "Y", zlab= "Z",
       xlim=c(-3000,3000), ylim=c(-3000,3000), zlim=c(-3000,3000))
```

### Spherical projections

#### g-sphere

A _g-sphere_ is a method of visualising tri-axial *acceleration* data. This involves centering the data and plotting it on a sphere.The function `triACC` allows the user to center this data (as well as calculating  pitch, roll and yaw from the data)

```r
# plot an g-phere
calibration = triACC(dta = PAM_data$magnetic)
pamSPHERE(x = calibration$centered_accx,
          y = calibration$centered_accy,
          z = calibration$centered_accz,
          ptcol = "royalblue",
          ptsize = 0.03,
          linecolor ="orange",
          spherecolor="orange",
          arrows=TRUE)


```

#### m-sphere

An _m-sphere_ is a method of visualising tri-axial *magnetometer* data. This involves centering the data and correcting the data, before  plotting it on a sphere.The function `triMAG` calibrates the data. This provides the animal's bearing.

```r
# plot a m-phere
calibration = triMAG(dta = PAM_data$magnetic)
pamSPHERE(x = calibration$calib_magx,
          y = calibration$calib_magy,
          z = calibration$calib_magz,
          ptcol = "orange",
          ptsize = 0.03,
          linecolor ="royalblue",
          spherecolor="royalblue",
          arrows=TRUE,
          cex=2)
```


## Classifying migratory flapping flight in Hoopoes  (an example classification)

## Performing the classification

Because flapping is widespread in birds, **PAMLr** integrates a pre-defined function `classifyFLAP()` to classify this behaviour.

This function assumes that if the bird has displayed **high activity** for x number of consecutive minutes, then it is flapping. It is therefore important to think about what constitutes high activity and how long this period should be. At the moment, the function uses **k-means clustering** to identify the threshold between high and low activity. Using `toPLOT = TRUE` then allows the user to see where that threshold was drawn. The period of high activity is set by default to `period = 3`. This is because activity is recorded (on this logger) every 5 minutes and we assume that after an hour of high activity, the bird must be flapping. 

Thus "high activity duration" / "data resolution" = "period" and 60 minutes / 5 minutes = period of 12.

```r
# Classify behaviour
behaviour = classifyFLAP(dta = PAM_data$acceleration, period = 12)
str(behaviour)
```
<img align="center" src="https://github.com/KiranLDA/PAMLrManual/blob/master/_bookdown_files/PAMLrManual_files/figure-html/unnamed-chunk-50-1.png">


This classification therefore provides different pieces of infomration.

* **timetable** shows when a migratory flapping flight started and stopped, and how long it lasted (in hours)
* **classification** is the output from the classification. In this case, each cluster/classs/state is represented by numbers between one 1 and 4. To find out what behaviour each of these numbers represent, we can refer to **low_movement**, **high_movement**, **migration** and **no_movement** 
* **threshold** represents the threshold between high and low activity.

Using these information, it's therefore possible to plot the classification:


```r
#Plot behaviour
col=col=c("black","royalblue4","brown","gold")
index= 7300:8000
plot(PAM_data$acceleration$date[index],PAM_data$acceleration$act[index],
  type="l", xlab="Date", ylab="Activity")
points(PAM_data$acceleration$date[index],PAM_data$acceleration$act[index],
  col=col[behaviour$classification+1][index], 
  pch=16,)
legend( PAM_data$acceleration$date[index[1]],60 , 
        c("No activity", "Low activity", "High activity", "Migration" ) ,
        col = col[c(behaviour$no_movement, behaviour$low_movement,
                    behaviour$high_movement, behaviour$migration)+1],
        pch=20)
```
<img align="center" src="https://github.com/KiranLDA/PAMLrManual/blob/master/docs/08-flapping_files/figure-html/unnamed-chunk-5-1.png">

## Plot the classification as a sensor image

Another way of looking at a classification is to use a sensor image of the results and to plot it side by side with the raw data to see if the same patterns are being picked out. We can also add (for instance sunset and sunrise events)

```r
par(mfrow= c(1,3), # number of panels
    oma=c(0,2,0,6), # outer margin around all panels
    mar =  c(4,1,4,1)) # inner margin around individual fivure

sensorIMG(PAM_data$acceleration$date, ploty=FALSE,
          PAM_data$acceleration$act, main = "Activity",
          col=c("black",viridis::cividis(90)), cex=1.2, cex.main = 2)
legend("bottomright",cex=1.2,
   c("No Activity", "Low Activity", "High Activity" ) , fill = c("black","royalblue3", "orange"), xpd = NA)


sensorIMG(PAM_data$pressure$date, plotx=TRUE, ploty=FALSE, labely=FALSE,
          PAM_data$pressure$obs,  main="Pressure",
          col=c("black",viridis::cividis(90)), cex=1.2, cex.main = 2)
legend("bottomright",cex=1.2,
   c("Low Pressure", "High Pressure" ) , fill = c("royalblue3", "orange"), xpd = NA)


sensorIMG(PAM_data$acceleration$date, labely=FALSE,
          behaviour$classification, 
          main="Classification",
          col=col,
          cex=1.2, cex.main = 2)
legend("bottomright",cex=1.2,
  # grconvertX(1, "device"), grconvertY(1, "device"),
   c("Resting", "Active", "Flapping", "Migrating" ) , fill = col, xpd = NA)


```

<img align="center" src="https://github.com/KiranLDA/PAMLrManual/blob/master/docs/08-flapping_files/figure-html/unnamed-chunk-6-1.png">



## Authors

Kiran Dhanjal-Adams

## Relevant references

Dhanjal-Adams, K.L., Bauer, S., Emmenegger, T., Hahn, S., Lisovski, S., & Liechti, F. (2018) [Spatiotemporal Group Dynamics in a Long-Distance Migratory Bird](https://www.cell.com/current-biology/fulltext/S0960-9822(18)30845-5?_returnURL=https%3A%2F%2Flinkinghub.elsevier.com%2Fretrieve%2Fpii%2FS0960982218308455%3Fshowall%3Dtrue). Current Biology, 28, 2824–2830.e3. 

Liechti, F., Bauer, S., Dhanjal-Adams, K.L., Emmenegger, T., Zehtindjiev, P., & Hahn, S. (2018) [Miniaturized multi-sensor loggers provide new insight into year-round flight behaviour of small trans-Sahara avian migrants](https://movementecologyjournal.biomedcentral.com/articles/10.1186/s40462-018-0137-1). Movement Ecology, 6, 19. 

## License

This project is licensed under the GNU General Public License version 3 - see the [LICENSE](https://github.com/KiranLDA/PAMLr/blob/master/LICENSE) file for details
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/create_rolling_window.R
\name{create_rolling_window}
\alias{create_rolling_window}
\title{Rolling window}
\usage{
create_rolling_window(
  dta,
  interp = TRUE,
  resolution_out = 15,
  window = 120,
  units = "mins"
)
}
\arguments{
\item{dta}{PAM data to be used in the analysis}

\item{interp}{By default TRUE. Whether or not to interpolate NAs in dataset that the rollapply is used on.}

\item{resolution_out}{Temporal resolution of output dataset. By defaukt 15 mins. Must be in minutes unless the units are changed}

\item{window}{Window over which to apply the rolling window. By defaukt 120 mins. Equivalent to zoo::rollapply(,width = window,) Must be in minutes unless the units are changed}

\item{units}{By default"mins", but also supports "hours" and "secs"}
}
\value{
a dataframe of derived metrics including the median, standard deviation, minimum, maximum, cumulative difference and range over the period
}
\description{
This function merges all the data to a given output resolution and then progresses along the timeseries and creates summary statistics within a pre-defined time window. Can interpolate or not. Interpolations is only recommended if the analysis cannot handle NAs and if the window is smaller than the data (e.g. magnetic) with the worst temporal resolution
}
\examples{

data(swift)
PAM_data = swift

# crop the data to get rid of no good periods
start = as.POSIXct("2016-09-01","\%Y-\%m-\%d", tz="UTC")
end = as.POSIXct("2017-04-21","\%Y-\%m-\%d", tz="UTC")
PAM_data = create_crop(PAM_data, start, end)

TOclassify = create_rolling_window(dta = list(pressure = PAM_data$pressure,
                                acceleration = PAM_data$acceleration),
                     resolution_out = 60 ,
                     window = 24*60)

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/triACC.R
\name{triACC}
\alias{triACC}
\title{Tri-axial acceleration functions}
\usage{
triACC(dta)
}
\arguments{
\item{dta}{magentic data from PAM logger see hoopoe$magnetic for an example}
}
\value{
roll, pitch and yaw from acceleration data
}
\description{
This function calculates roll, pitch and yaw from tri-axial acceleration data. It also centeres the acceleration data (for plotting on a sphere).
}
\examples{
#data(swift)
#PAM_data = swift

#triACC(dta = PAM_data$magnetic)

}
\references{
Bidder, O.R., Walker, J.S., Jones, M.W., Holton, M.D., Urge, P., Scantlebury, D.M., Marks, N.J., Magowan, E.A., Maguire, I.E. and Wilson, R.P., 2015. Step by step: reconstruction of terrestrial animal movement paths by dead-reckoning. Movement ecology, 3(1), p.23.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rollPAM.R
\name{rollPAM}
\alias{rollPAM}
\title{Rolling window}
\usage{
rollPAM(dta, interp = TRUE, resolution_out = 15, window = 120, units = "mins")
}
\arguments{
\item{dta}{PAM data to be used in the analysis}

\item{interp}{By default TRUE. Whether or not to interpolate NAs in dataset that the rollapply is used on.}

\item{resolution_out}{Temporal resolution of output dataset. By defaukt 15 mins. Must be in minutes unless the units are changed}

\item{window}{Window over which to apply the rolling window. By defaukt 120 mins. Equivalent to zoo::rollapply(,width = window,) Must be in minutes unless the units are changed}

\item{units}{By default"mins", but also supports "hours" and "secs"}
}
\value{
a dataframe of derived metrics including the median, standard deviation, minimum, maximum, cumulative difference and range over the period
}
\description{
This function merges all the data to a given output resolution and then progresses along the timeseries and creates summary statistics within a pre-defined time window. Can interpolate or not. Interpolations is only recommended if the analysis cannot handle NAs and if the window is smaller than the data (e.g. magnetic) with the worst temporal resolution
}
\examples{

#data(swift)
#PAM_data = swift

## crop the data to get rid of no good periods
#start = as.POSIXct("2016-09-01","\%Y-\%m-\%d", tz="UTC")
#end = as.POSIXct("2017-04-21","\%Y-\%m-\%d", tz="UTC")
#PAM_data = cutPAM(PAM_data, start, end)

#TOclassify = rollPAM(dta = list(pressure = PAM_data$pressure,
#                                acceleration = PAM_data$acceleration),
#                     resolution_out = 60 ,
#                     window = 24*60)

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/reducePAM.R
\name{reducePAM}
\alias{reducePAM}
\title{Reduces the data to a specified temporal resolution}
\usage{
reducePAM(dta, varint, interp = TRUE, summary = "median")
}
\arguments{
\item{dta}{raw pam data see \verb{data(bee_eater} for example}

\item{varint}{the variable of interest.Sipprots c("pressure","light","acceleration","temperature","magnetic"). For instance if all other variables should be summarised at the same temporal resolution as this specified variable.}

\item{interp}{Default TRUE. Whether or not to replace NAs with interpolated values}

\item{summary}{Can be "sum", "median" or "none". What type of summary variable to give when condensing data -}
}
\value{
reduced/summarised and interpolated dataset
}
\description{
This function takesthe typical PAM_data input, which is a nested list of different sensor data, all formatted at different time resolutions. By specifying a particular sensor, all data are formatted at the same temporal resolution as this data. By default, the median of data with smaller temporal resolution are kept, and data are interpolated.
}
\examples{
#data(bee_eater)
#PAM_data = bee_eater

#start = as.POSIXct("2015-08-01","\%Y-\%m-\%d", tz="UTC")
#end = as.POSIXct("2016-06-21","\%Y-\%m-\%d", tz="UTC")
#PAM_data = cutPAM(PAM_data, start, end)

#reduced_dta = reducePAM(PAM_data , "pressure", interp = FALSE)
#head(reduced_dta)

#reduced_dta = reducePAM(PAM_data , "act", interp = FALSE)
#head(reduced_dta)

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_sensorimage.R
\name{plot_sensorimage}
\alias{plot_sensorimage}
\title{sensor image}
\usage{
plot_sensorimage(
  date,
  sensor_data,
  tz = "UTC",
  plotx = TRUE,
  ploty = TRUE,
  labelx = TRUE,
  labely = TRUE,
  offset = 0,
  dt = NA,
  xlab = "Hour",
  ylab = "Date",
  cex = 2,
  na.col = "white",
  col = c("black", viridis::magma(90)),
  ...
)
}
\arguments{
\item{date}{Date data in POSIXct format, most commonly \code{PAM_data$acceleration$date}}

\item{sensor_data}{sensor data, for example look at \code{PAM_data$acceleration$act}}

\item{tz}{Time zone for POSIXct, default set to "UTC"}

\item{plotx}{wherether or not to plot the x axis ticks + labels (for instance when compiling multifigures)}

\item{ploty}{wherether or not to plot the y axis ticks + labels (for instance when compiling multifigures)}

\item{labelx}{wherether or not to write the name of the x axis (for instance when compiling multifigures)}

\item{labely}{wherether or not to write the name of the y axis (for instance when compiling multifigures)}

\item{offset}{This parameter determines where the center of the graph is. When \code{offset = 0}, then midday is at the center of the graph. when \code{offset=12} midnight`}

\item{dt}{the time interval to which the data are resampled (secs). Default is \code{NA}}

\item{xlab}{label for x-axis (as a character string)}

\item{ylab}{label for y axis (as a character string)}

\item{cex}{size of labels}

\item{na.col}{colour given to NA values, default is "white"}

\item{col}{Colour scheme of plot. Default \code{col = c("black",viridis::magma(90))}}

\item{...}{Any additional parameters used by graphics::image}
}
\value{
an image of the sensor data, for instance with activity it would produce an actogram
}
\description{
This function plots sensor data as an image. An actogram for instance is a type of sensor image.
}
\examples{
#specify the data location
data(hoopoe)
start = as.POSIXct("2016-07-01","\%Y-\%m-\%d", tz="UTC")
end = as.POSIXct("2017-06-01","\%Y-\%m-\%d", tz="UTC")
PAM_data = create_crop(hoopoe,start,end)

# Create plots with 3 together (mfrow)
par( mfrow= c(1,3), oma=c(0,2,0,6))

par(mar =  c(4,2,4,2))
plot_sensorimage(PAM_data$acceleration$date, ploty=FALSE,
          PAM_data$acceleration$act, main = "Activity",
          col=c("black",viridis::cividis(90)), cex=1.2, cex.main = 2)

par(mar =  c(4,2,4,2))
plot_sensorimage(PAM_data$pressure$date, plotx=TRUE, ploty=FALSE, labely=FALSE,
          PAM_data$pressure$obs,  main="Pressure",
          col=c("black",viridis::cividis(90)), cex=1.2, cex.main = 2)

par(mar =  c(4,2,4,2))
plot_sensorimage(PAM_data$temperature$date, labely=FALSE,
          PAM_data$temperature$obs,  main="Temperature",
          col=c("black",viridis::cividis(90)), cex=1.2, cex.main = 2)

######################################################
# Look at a classification output
######################################################

# Classification
classification  =  classify_flap(dta = PAM_data$acceleration, period = 10, to_plot=FALSE)

par( mfrow= c(1,3), oma=c(0,2,0,6),mar =  c(4,2,4,2))

plot_sensorimage(PAM_data$pressure$date, c(0,abs(diff(PAM_data$pressure$obs))),
          main="Pressure  difference",
          ploty=FALSE,
          col=c("black",viridis::cividis(90)), cex=1.2, cex.main = 2)

plot_sensorimage(PAM_data$acceleration$date, PAM_data$acceleration$act,  main="Activity",
          ploty=FALSE, labely=FALSE,
          col=c(viridis::cividis(90)), cex=1.2, cex.main = 2)

plot_sensorimage(PAM_data$acceleration$date,
          ifelse(classification$classification == classification$migration, 1,2),
          main="Migration Classification",
          labely=FALSE,
          na.col="white",
          col = c("orange","black"),
          cex=1.2, cex.main = 2)


twilights <- GeoLight::twilightCalc(PAM_data$light$date,
                                    PAM_data$light$obs,
                                    LightThreshold = 2,
                                    ask = FALSE)

plot_sensorimage_twilight(twilights$tFirst, offset=0,
       col= ifelse(twilights$type == 1,
                   "goldenrod","cornflowerblue"),
       pch=16, cex=0.5)


}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/create_import.R
\name{create_import}
\alias{create_import}
\title{Import PAM data}
\usage{
create_import(
  pathname = pathname,
  measurements = c(".pressure", ".glf", ".acceleration", ".temperature", ".magnetic")
)
}
\arguments{
\item{pathname}{path where files are stored}

\item{measurements}{a series of measurements logged by the PAM logger which are to be imported. Currently supports these file extentions: ".pressure", ".glf", ".gle",".acceleration", ".temperature", "AirTemperature", ".BodyTemperature" and ".magnetic"}
}
\value{
a list of measurements for the one individual
}
\description{
Imports and formats many datasets into one big nested list containing all the data from the different sensors. A subset of sensors can be selected using \code{measurements}.
}
\examples{
#pathname = "your/path/here"
#measurements = c(".pressure", ".glf")
#PAM_data = create_import(pathname, measurements)
#str(PAM_data)
#plot(PAM_data$light$date[3000:5000], PAM_data$light$obs[3000:5000],
#type="l", xlab="Date", ylab="Light Intensity")

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/create_summary_statistics.R
\name{create_summary_statistics}
\alias{create_summary_statistics}
\title{Prepare data for analysis by deriving summary statistics}
\usage{
create_summary_statistics(
  dta,
  availavariable = c("pressure", "light", "acceleration", "magnetic", "temperature"),
  Pdiff_thld = 2,
  light_thld = 1,
  method = "pressure",
  twl,
  interp = FALSE,
  tz = "UTC"
)
}
\arguments{
\item{dta}{PAM data to be used in the analysis e.g. str(hoopoe)}

\item{availavariable}{Variables to be used to derive metrics for classification. must have "pressure", but ideally \code{availavariable = c("pressure", "light", "acceleration")} if any of these are incomplete, do not use them}

\item{Pdiff_thld}{Pressure threshold. Only used when method="pressure".  This if pressure changes more than e.g. 2hpa over 30 minutes, then the bird is flying}

\item{light_thld}{Light threshold. Only used when method="darkness". This is the the light threshold for finding darkness, should be the same as for GeoLight::twilightCalc}

\item{method}{The type of event that is being classified. Can be "flap", "endurance", "rest", "pressure","light" or "darkness".If method = "pressure" then  it find periods where pressure has changed more than a certain threshold. If method = "flap", then the algorithm looks for sustained periods of high activity. If method = "endurance" it looks for sustained activity (low or high). If method = "rest+ then it looks for sustained periods of no activity. If method = "light" if looks for periods of sustained sunlight. If method = "darkness" if looks for periods of darkness}

\item{twl}{twilight estimates formatted according to twilightCalc in GeoLight package}

\item{interp}{whether or not to interpolate the magnetic data. If FALSE, then NAs are left in the dataset}

\item{tz}{Timeuzone. default is "UTC"}
}
\value{
a dataframe of derives metrhics based on pressure:

date : Date (without time)

start : Start time and date of the event, POSIXct format

end : Time and date that the event finished, POSIXct format

duration : How long it lasted (in hours)

total_daily_duration : The total duration of all the events that occured that day (in hours)

total_daily_event_number : The total number of events which occured that day

cum_pressure_change : The cumulative change in atmospheric pressure during that event (in hectopascals)

cum_altitude_change : The cumulative change in altitude during that event (in metres)

cum_altitude_up : The cumulative number of metres that the bird went upwards during that event

total_daily_P_change : The cumulative change in atmospheric pressure for all the events for that date (in hectopascals)

P_dep_arr : The difference between atmospheric pressure at the start of the event, and at the end (in hectopascals)

pressure_range : The total range of the atmospheric pressure during that event (maximum minus miniimum - in hectopascals)

altitude_range : The total altitude range during that event (maximum minus miniimum - in metres)

mean_night_P : The mean pressure during the night before the event took place (in hectopascals)

sd_night_P : The standard deviation of pressure the night before the event took place (in hectopascals)

mean_nextnight_P : The mean pressure the night after the event took place (in hectopascals)

sd_nextnight_P : The standard deviation of pressure the night after the event took place (in hectopascals)

night_P_diff : The difference between the mean pressures of the night before and the night after the event took place (in hectopascals)

median_activity : The median activity during that event

sum_activity : The sum of the activity during that event

prop_resting : The propotion of time during that event where activity = 0

prop_active : The propotion of time during that event where activity > 0

mean_night_act : The mean activity during the night before the event took place

sd_night_act : The standard deviation of activity the night before the event took place

sum_night_act : The summed activity during the night before the event took place

mean_nextnight_act :The mean activity the night after the event took place

sd_nextnight_act : The standard deviation of activity the night after the event took place

sum_nextnight_act : The summed activity the night after the event took place

night_act_diff : The difference between the mean activity of the night before and the night after the event took place

median_pitch : The median pitch during that event

sd_pitch : The standard deviation of pitch during that event

median_light : The median light recordings during that event

nightime : Whether or not it was night during the majority of the event (1= night, 0 = day)

median_gX : Median raw acceledation on the x axis during the event

median_gY : Median raw acceledation on the y axis during the event

median_gZ : Median raw acceledation on the z axis during the event

median_mX : Median raw magnetic field on the x axis during the event

median_mY : Median raw magnetic field on the y axis during the event

median_mZ : Median raw magnetic field on the z axis

median_temp : Median temperature during the event (in celsius)

sd_temp : Standard deviation of temperature during the event (in celsius)

cum_temp_change : Cumulative absolute difference in temperature during the event (in celsius)
}
\description{
This function summarises the data based on different patterns such as sustained acitivity, sustained pressure changes, etc... and extracts these from the data as a timetable. It then creates summary statics for each of these periods or events, such as cumulative altitude change, mean pitch etc...
}
\examples{
data(hoopoe)
PAM_data=hoopoe
twl = GeoLight::twilightCalc(PAM_data$light$date, PAM_data$light$obs,
                             LightThreshold = 2, ask = FALSE)

to_classify = create_summary_statistics(PAM_data,
                     method= "flap",
                     twl = twl)

str(to_classify)

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/compareCLASS.R
\name{compareCLASS}
\alias{compareCLASS}
\title{Combine classification results}
\usage{
compareCLASS(date, classifications)
}
\arguments{
\item{date}{datetime in POSIXCT format. See hoopoe$activity$date for example}

\item{classifications}{a dataframe containing the results from different classifications in each column. The classifications must correspond to the same datetimes.}
}
\value{
a dataframe containing date

the raw classifications

each state with the number of times it was used in a classification

whether or not all classes provided the same state
}
\description{
this function takes multiple classifications of the same data and compares the agreement between each, on a timepoint level.
}
\examples{

## Import data
#data(hoopoe)
#start = as.POSIXct("2016-07-01","\%Y-\%m-\%d", tz="UTC")
#end = as.POSIXct("2017-06-01","\%Y-\%m-\%d", tz="UTC")
#PAM_data= cutPAM(hoopoe,start,end)

## perform one classification using classifyFLAP
#classification = classifyFLAP(dta = PAM_data$acceleration, period = 12)
## Put the classification in the same resolution as pressure
#class1 = classification2PAM(from = classification$timetable$start,
#                            to = classification$timetable$end,
#                            # because the timetable only contains migration periods
#                            classification = rep_len(1,length(classification$timetable$end)),
#                            addTO = PAM_data$pressure)
## Convert to categories
#class1 = ifelse(class1 == classification$migration, "Migration", "Other")


## Perform another classification using pressure difference
#class2 = c(0,ifelse(abs(diff(PAM_data$pressure$obs))>2, "Migration", "Other"))

## both classes have been converted to the same time intervals as pressure,
#date = PAM_data$pressure$date

## Combine the classifications into a dataframe
#classifications = data.frame(flap = class1, # flapping classification
#                            Pdiff = class2) # pressure difference classification

#class_comparison = compareCLASS(date=date,
#                                classifications=classifications)

#plot(PAM_data$pressure$date,
#     PAM_data$pressure$obs,
#     type = "l",
#     xlab = "Date",
#     ylab = "Pressure (hPa)",
#     col = "royalblue3")

#points(PAM_data$pressure$date,
#       PAM_data$pressure$obs,
#      cex = class_comparison$Migration / 2,
#       col ="orange",
#       pch = 16)

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calculate_triaxial_acceleration.R
\name{calculate_triaxial_acceleration}
\alias{calculate_triaxial_acceleration}
\title{Tri-axial acceleration functions}
\usage{
calculate_triaxial_acceleration(dta)
}
\arguments{
\item{dta}{magentic data from PAM logger see hoopoe$magnetic for an example}
}
\value{
roll, pitch and yaw from acceleration data
}
\description{
This function calculates roll, pitch and yaw from tri-axial acceleration data. It also centeres the acceleration data (for plotting on a sphere).
}
\examples{
data(swift)
PAM_data = swift

calculate_triaxial_acceleration(dta = PAM_data$magnetic)

}
\references{
Bidder, O.R., Walker, J.S., Jones, M.W., Holton, M.D., Urge, P., Scantlebury, D.M., Marks, N.J., Magowan, E.A., Maguire, I.E. and Wilson, R.P., 2015. Step by step: reconstruction of terrestrial animal movement paths by dead-reckoning. Movement ecology, 3(1), p.23.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/altitudeCALC.R
\name{altitudeCALC}
\alias{altitudeCALC}
\title{Calculate altitude from pressure}
\usage{
altitudeCALC(P, T0 = 288.15, P0 = 1013.25)
}
\arguments{
\item{P}{on-bird pressure in hectopascals}

\item{T0}{temperature at sea level in kelvin (i.e. celcius + 273.15 OR (Farenheit − 32) * 5/9 + 273.15)}

\item{P0}{pressure at sea level in hectopascals}
}
\value{
altitude in metres, default is calculated according to International Standard Atmosphere model (International Organization for Standardization 1975: ISO 2533:1975)
}
\description{
This function calculated altitude from atmospheric pressure
}
\references{
Atmosphere, S., 1975. International organization for standardization. ISO, 2533, p.1975.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/classify_changepoint.R
\name{classify_changepoint}
\alias{classify_changepoint}
\title{Changepoint analysis}
\usage{
classify_changepoint(
  dta,
  cpt.method = "meanvar",
  method = "PELT",
  penalty = "Manual",
  pen.value = "100*log(n)",
  ...
)
}
\arguments{
\item{dta}{data to be classified (vector)}

\item{cpt.method}{method for classifying data, currently support "mean","variance", "meanvar", use ?changepoint::cpt.mean ?changepoint::cpt.var or ?changepoint::cpt.meanvar for additional parameters}

\item{method}{changepoint method used in package changepoint. see changepoint::cpt.mean, changepoint::cpt.var and changepoint::cpt.meanvar for details}

\item{penalty}{See changepoint  package for details. Choice of "None", "SIC", "BIC", "MBIC", AIC", "Hannan-Quinn", "Asymptotic", "Manual" and "CROPS" penalties. If Manual is specified, the manual penalty is contained in the pen.value parameter. If Asymptotic is specified, the theoretical type I error is contained in the pen.value parameter. If CROPS is specified, the penalty range is contained in the pen.value parameter; note this is a vector of length 2 which contains the minimum and maximum penalty value. Note CROPS can only be used if the method is "PELT". The predefined penalties listed DO count the changepoint as a parameter, postfix a 0 e.g."SIC0" to NOT count the changepoint as a parameter.}

\item{pen.value}{See changepoint package for details. Tthe theoretical type I error e.g.0.05 when using the Asymptotic penalty. A vector of length 2 (min,max) if using the CROPS penalty. The value of the penalty when using the Manual penalty option - this can be a numeric value or text giving the formula to use. Available variables are, n=length of original data, null=null likelihood, alt=alternative likelihood, tau=proposed changepoint, diffparam=difference in number of alternatve and null parameters}

\item{...}{any additional inputs for changepoint::cpt.mean,  changepoint::cpt.var or changepoint::cpt.meanvar}
}
\value{
Changepoints in the data.
}
\description{
This function implements a changepoint analysis to find changes in variance, mean or both. By default it is set to look for changes in the pressure data to identify migration periods, however it can be altered to any form of changepoint analysis in the package changepoint.
}
\examples{

# Import and crop PAM data
data("swift")
start = as.POSIXct("2016-09-01","\%Y-\%m-\%d", tz="UTC")
end = as.POSIXct("2017-04-15","\%Y-\%m-\%d", tz="UTC")
PAM_data =  create_crop(swift, start, end)

# data(bee_eater)
# start = as.POSIXct("2015-08-01","\%Y-\%m-\%d", tz="UTC")
# end = as.POSIXct("2016-06-01","\%Y-\%m-\%d", tz="UTC")
# PAM_data = create_crop(bee_eater, start, end)

changepoints  = classify_changepoint(PAM_data$pressure$obs)

# plot the timeseries with the changepoint
par(mfrow=c(2,1))
plot(PAM_data$pressure$obs, type="l")
abline(v=changepoints$changepoints, col="red",lwd=2)

# plot using changepoint package output
changepoint::plot(changepoints$output, cpt.width=3)

}
\references{
Killick, R. and Eckley, I., 2014. changepoint: An R package for changepoint analysis. Journal of statistical software, 58(3), pp.1-19.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_interactive_timeseries.R
\name{plot_interactive_timeseries}
\alias{plot_interactive_timeseries}
\title{Plot PAM data as an interactive timeseries}
\usage{
plot_interactive_timeseries(
  dta,
  from = dta$light$date[1],
  to = dta$light$date[length(dta$light$date)],
  to_plot = names(dta)
)
}
\arguments{
\item{dta}{PAM data to be plotted}

\item{from}{date that plotting starts}

\item{to}{date that plotting ends}

\item{to_plot}{names of the variables to plot. For now this includes \code{light}, \code{pressure}, \code{acceleration} and \code{temperature}}
}
\value{
a plot of all the measurements
}
\description{
This opens a java dygraph application which allows the user to zoom in and out. In Rstudio it will open in the viewer pane and in base R in an html readers. Note that this can be a bit slow
}
\examples{
#load dummy data
data(hoopoe)
PAM_data=hoopoe

# This bit is for Rstudio users to prevent html from opening in Viewer pane and crashing
# It opens in web browser instead
backup_options <- options()
options(viewer=NULL)

# Plot interactive graphics
plot_interactive_timeseries(dta = PAM_data)

# restore Rstudio settings from before plot
options(backup_options)


}
\references{
Vanderkam, D., Allaire, J., Owen, J., Gromer, D., Shevtsov, P. and Thieurmel, B., dygraphs: Interface to Dygraphs Interactive Time Series Charting Library, 2015. URL http://CRAN. R-project. org/package= dygraphs. R package version 0.4, 5, p.7.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/convert_diveMove.R
\name{convert_diveMove}
\alias{convert_diveMove}
\title{convert diveMove input into pamlr input}
\usage{
convert_diveMove(
  input,
  id = 0,
  measurements = c("light", "depth", "temperature")
)
}
\arguments{
\item{input}{diveMove data}

\item{id}{individual identifyer for tracked animal e.g. leg tag number}

\item{measurements}{a series of measurements logged by the dtr logger which are to be imported. Currently supports: "light", "depth","temperature"}
}
\value{
a list of measurements for the one individual
}
\description{
Converts  diveMove input into pamlr input
}
\examples{
#install.packages("diveMove")
#library(diveMove)
#zz <- system.file(file.path("data", "dives.csv"),
#                  package="diveMove", mustWork=TRUE)
#(sealX <- readTDR(zz, speed=TRUE, sep=";", na.strings="", as.is=TRUE))

#PAM_data = convert_diveMove(sealX, id="sealX")

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pamPREP.R
\name{pamPREP}
\alias{pamPREP}
\title{Prepare data for analysis by deriving summary statistics}
\usage{
pamPREP(
  dta,
  availavariable = c("pressure", "light", "acceleration", "magnetic", "temperature"),
  Pdiff_thld = 2,
  light_thld = 1,
  method = "pressure",
  twl,
  interp = FALSE,
  tz = "UTC"
)
}
\arguments{
\item{dta}{PAM data to be used in the analysis e.g. str(hoopoe)}

\item{availavariable}{Variables to be used to derive metrics for classification. must have "pressure", but ideally \code{availavariable = c("pressure", "light", "acceleration")} if any of these are incomplete, do not use them}

\item{Pdiff_thld}{Pressure threshold. Only used when method="pressure".  This if pressure changes more than e.g. 2hpa over 30 minutes, then the bird is flying}

\item{light_thld}{Light threshold. Only used when method="darkness". This is the the light threshold for finding darkness, should be the same as for GeoLight::twilightCalc}

\item{method}{The type of event that is being classified. Can be "flap", "endurance", "rest", "pressure","light" or "darkness".If method = "pressure" then  it find periods where pressure has changed more than a certain threshold. If method = "flap", then the algorithm looks for sustained periods of high activity. If method = "endurance" it looks for sustained activity (low or high). If method = "rest+ then it looks for sustained periods of no activity. If method = "light" if looks for periods of sustained sunlight. If method = "darkness" if looks for periods of darkness}

\item{twl}{twilight estimates formatted according to twilightCalc in GeoLight package}

\item{interp}{whether or not to interpolate the magnetic data. If FALSE, then NAs are left in the dataset}

\item{tz}{Timeuzone. default is "UTC"}
}
\value{
a dataframe of derives metrhics based on pressure:

date : Date (without time)

start : Start time and date of the event, POSIXct format

end : Time and date that the event finished, POSIXct format

duration : How long it lasted (in hours)

total_daily_duration : The total duration of all the events that occured that day (in hours)

total_daily_event_number : The total number of events which occured that day

cum_pressure_change : The cumulative change in atmospheric pressure during that event (in hectopascals)

cum_altitude_change : The cumulative change in altitude during that event (in metres)

cum_altitude_up : The cumulative number of metres that the bird went upwards during that event

total_daily_P_change : The cumulative change in atmospheric pressure for all the events for that date (in hectopascals)

P_dep_arr : The difference between atmospheric pressure at the start of the event, and at the end (in hectopascals)

pressure_range : The total range of the atmospheric pressure during that event (maximum minus miniimum - in hectopascals)

altitude_range : The total altitude range during that event (maximum minus miniimum - in metres)

mean_night_P : The mean pressure during the night before the event took place (in hectopascals)

sd_night_P : The standard deviation of pressure the night before the event took place (in hectopascals)

mean_nextnight_P : The mean pressure the night after the event took place (in hectopascals)

sd_nextnight_P : The standard deviation of pressure the night after the event took place (in hectopascals)

night_P_diff : The difference between the mean pressures of the night before and the night after the event took place (in hectopascals)

median_activity : The median activity during that event

sum_activity : The sum of the activity during that event

prop_resting : The propotion of time during that event where activity = 0

prop_active : The propotion of time during that event where activity > 0

mean_night_act : The mean activity during the night before the event took place

sd_night_act : The standard deviation of activity the night before the event took place

sum_night_act : The summed activity during the night before the event took place

mean_nextnight_act :The mean activity the night after the event took place

sd_nextnight_act : The standard deviation of activity the night after the event took place

sum_nextnight_act : The summed activity the night after the event took place

night_act_diff : The difference between the mean activity of the night before and the night after the event took place

median_pitch : The median pitch during that event

sd_pitch : The standard deviation of pitch during that event

median_light : The median light recordings during that event

nightime : Whether or not it was night during the majority of the event (1= night, 0 = day)

median_gX : Median raw acceledation on the x axis during the event

median_gY : Median raw acceledation on the y axis during the event

median_gZ : Median raw acceledation on the z axis during the event

median_mX : Median raw magnetic field on the x axis during the event

median_mY : Median raw magnetic field on the y axis during the event

median_mZ : Median raw magnetic field on the z axis

median_temp : Median temperature during the event (in celsius)

sd_temp : Standard deviation of temperature during the event (in celsius)

cum_temp_change : Cumulative absolute difference in temperature during the event (in celsius)
}
\description{
This function summarises the data based on different patterns such as sustained acitivity, sustained pressure changes, etc... and extracts these from the data as a timetable. It then creates summary statics for each of these periods or events, such as cumulative altitude change, mean pitch etc...
}
\examples{
#data(hoopoe)
#PAM_data=hoopoe
#twl = GeoLight::twilightCalc(PAM_data$light$date, PAM_data$light$obs,
#                             LightThreshold = 2, ask = FALSE)

#TOclassify = pamPREP(PAM_data,
#                     method= "flap",
#                     twl = twl)

#str(TOclassify)

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_interactive_sphere.R
\name{plot_interactive_sphere}
\alias{plot_interactive_sphere}
\title{Plot 3d sphere}
\usage{
plot_interactive_sphere(
  x,
  y,
  z,
  spherecolor = "white",
  linecolor = "black",
  linewidth = 2,
  ptcol = "black",
  ptsize = 0.02,
  arrows = TRUE,
  cex = 1.5,
  ...
)
}
\arguments{
\item{x}{data on x axis}

\item{y}{data on y axis}

\item{z}{data on z axis}

\item{spherecolor}{color of the sphere}

\item{linecolor}{color of the lines on the sphere}

\item{linewidth}{width of the lines on the sphere}

\item{ptcol}{color of points}

\item{ptsize}{size of the points plotted on the sphere}

\item{arrows}{default TRUE - whether or not to draw the arrows}

\item{cex}{size of font. If cex = 2 then fond is two times bigger.}

\item{...}{additional input for rgl::spheres3d for plotting points on the sphere using}
}
\value{
a 3d sphere plot
}
\description{
This function takes a spherical projection of acceleration or magentic data and plots it on a sphere which can be turned and zoomed into.
}
\examples{
data("swift")
start = as.POSIXct("2016-09-01","\%Y-\%m-\%d", tz="UTC")
end = as.POSIXct("2017-04-15","\%Y-\%m-\%d", tz="UTC")
swift = create_crop(swift, start, end)
PAM_data = swift

# plot an m-phere
calibration = calculate_triaxial_magnetic(dta = PAM_data$magnetic)
plot_interactive_sphere(x = calibration$calib_magx,
          y = calibration$calib_magy,
          z = calibration$calib_magz,
          ptcol = "goldenrod",
          ptsize = 0.03,
          linecolor ="black",
          spherecolor="royalblue4",
          arrows=TRUE,
          cex=2)

# plot an g-phere
calibration = calculate_triaxial_acceleration(dta = PAM_data$magnetic)
plot_interactive_sphere(x = calibration$centered_accx,
          y = calibration$centered_accy,
          z = calibration$centered_accz,
          ptcol = "royalblue4",
          ptsize = 0.03,
          linecolor ="black",
          spherecolor="goldenrod",
          arrows=TRUE)

}
\references{
Adler, D., Nenadic, O. and Zucchini, W., 2003, March. Rgl: A r-library for 3d visualization with opengl. In Proceedings of the 35th Symposium of the Interface: Computing Science and Statistics, Salt Lake City (Vol. 35).
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calculate_altitude.R
\name{calculate_altitude}
\alias{calculate_altitude}
\title{Calculate altitude from pressure}
\usage{
calculate_altitude(P, T0 = 288.15, P0 = 1013.25)
}
\arguments{
\item{P}{on-bird pressure in hectopascals}

\item{T0}{temperature at sea level in kelvin (i.e. celcius + 273.15 OR (Farenheit − 32) * 5/9 + 273.15)}

\item{P0}{pressure at sea level in hectopascals}
}
\value{
altitude in metres calculated according to International Standard Atmosphere model (International Organization for Standardization 1975: ISO 2533:1975)
}
\description{
This function calculates altitude from atmospheric pressure
}
\examples{
altitude = calculate_altitude(P = hoopoe$pressure$obs)
plot(hoopoe$pressure$date,
     altitude,
     type="o",
     pch=16,
     xlab="Date",
     ylab="Altitude (m)")

}
\references{
Atmosphere, S., 1975. International organization for standardization. ISO, 2533, p.1975.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/classifyPAM.R
\name{classifyPAM}
\alias{classifyPAM}
\title{Derive classification data for soaring birds}
\usage{
classifyPAM(
  dta,
  method = c("kmeans", "hmm", "embc", "agnes", "diana"),
  states = 2,
  family = stats::gaussian(),
  ...
)
}
\arguments{
\item{dta}{data to be classified, can be anything}

\item{method}{method for classifying data, currently support "hmm" for hidden markov model and "kmeans" for kmeans clustering}

\item{states}{number of states to classify the data into}

\item{family}{By default \code{gaussian()}. Parameter only used in the hmm method. Must be the same as the parameters used by \code{depmixS4::depmix}}

\item{...}{any additional inputs for depmixs4::depmix, stats::kmeans, cluster::diana, cluster::diana  or EMbC::embc functions, depending on method selected}
}
\value{
the data's classification based on the chosen algorithm
}
\description{
This function takes data, typically a dataframe of events, which are classified into different behavioural states using different algorithms.
}
\examples{

#######################################################
#  data prep
#######################################################

#data(bee_eater)
#start = as.POSIXct("2015-07-01","\%Y-\%m-\%d", tz="UTC")
#end = as.POSIXct("2016-06-01","\%Y-\%m-\%d", tz="UTC")
#PAM_data = cutPAM(bee_eater, start, end)

# twl = GeoLight::twilightCalc(PAM_data$light$date, PAM_data$light$obs,
#LightThreshold = 2, ask = FALSE)
#availavariable = c("pressure", "light", "acceleration")

#TOclassify = pamPREP(PAM_data,
#                      method= "flap",
#                      twl = twl)

#TOclassify = TOclassify[complete.cases(TOclassify),]

#######################################################
# k-means example
#######################################################

#classification = classifyPAM(TOclassify[,c("cum_altitude_change", "night_P_diff" )],
#                             states=2, "kmeans")$cluster
#pressure_classification = classification2PAM(from = TOclassify$start,
#to =TOclassify$end,
#classification = classification,
#addTO = PAM_data$pressure)

#plot(PAM_data$pressure$date, PAM_data$pressure$obs,
#     type="l")
#points(PAM_data$pressure$date, PAM_data$pressure$obs,
#       col= pressure_classification+1,
#       pch=16)

#######################################################
# HMM example
#######################################################

#classification = classifyPAM(TOclassify[,c("cum_altitude_change", "night_P_diff" )]
#                             #TOclassify$night_P_diff +
#                             #TOclassify$cum_altitude_change
#                             ,
#                             states=2, "hmm")$cluster
#pressure_classification = classification2PAM(from = TOclassify$start,
#to =TOclassify$end,
#classification = classification,
#addTO = PAM_data$pressure)

#plot(PAM_data$pressure$date, PAM_data$pressure$obs,
#     type="l")
#points(PAM_data$pressure$date, PAM_data$pressure$obs,
#       col= pressure_classification+1,
#       pch=16)


#######################################################
# EMBC example
#######################################################

#classification = classifyPAM(TOclassify[,c("cum_altitude_change", "night_P_diff" )],
#                             "embc")

#pressure_classification = classification2PAM(from = TOclassify$start,
#                                             to =TOclassify$end,
#                                             classification = classification$cluster,
#                                             addTO = PAM_data$pressure)

#plot(PAM_data$pressure$date, PAM_data$pressure$obs,
#     type="l")
#points(PAM_data$pressure$date, PAM_data$pressure$obs,
#       col= classification$output@C[pressure_classification+1],
#       pch=16)


#######################################################
# agnes example
#######################################################

#classification = classifyPAM(TOclassify[,c("cum_altitude_change", "night_P_diff" )],
#                             states = 2,
#                            "agnes")

#plot(classification$output, main="agnes", which.plot = 2)

#pressure_classification = classification2PAM(from = TOclassify$start,
#                                             to =TOclassify$end,
#                                             classification = classification$cluster,
#                                             addTO = PAM_data$pressure)

#plot(PAM_data$pressure$date, PAM_data$pressure$obs,
#     type="l")
#points(PAM_data$pressure$date, PAM_data$pressure$obs,
#       col= pressure_classification+1,
#       pch=16)

#######################################################
# diana example
#######################################################

#classification = classifyPAM(TOclassify[,c("cum_altitude_change", "night_P_diff" )],
#                             states = 2,
#                             "diana")

#plot(classification$output, which.plot = 2, main="diana")

#pressure_classification = classification2PAM(from = TOclassify$start,
#                                             to =TOclassify$end,
#                                             classification = classification$cluster,
#                                             addTO = PAM_data$pressure)

#plot(PAM_data$pressure$date, PAM_data$pressure$obs,
#     type="l")
#points(PAM_data$pressure$date, PAM_data$pressure$obs,
#       col= pressure_classification+1,
#       pch=16)



}
\references{
Forgy, E. W. (1965). Cluster analysis of multivariate data: efficiency vs interpretability of classifications. Biometrics, 21, 768–769.

Hartigan, J. A. and Wong, M. A. (1979). Algorithm AS 136: A K-means clustering algorithm. Applied Statistics, 28, 100–108. doi: 10.2307/2346830.

Lloyd, S. P. (1957, 1982). Least squares quantization in PCM. Technical Note, Bell Laboratories. Published in 1982 in IEEE Transactions on Information Theory, 28, 128–137.

MacQueen, J. (1967). Some methods for classification and analysis of multivariate observations. In Proceedings of the Fifth Berkeley Symposium on Mathematical Statistics and Probability, eds L. M. Le Cam & J. Neyman, 1, pp. 281–297. Berkeley, CA: University of California Press.

Ingmar Visser and Maarten Speekenbrink (2010). depmixS4: An R Package for Hidden Markov Models. Journal of Statistical Software, 36(7), p. 1-21.

Lawrence R. Rabiner (1989). A tutorial on hidden Markov models and selected applications in speech recognition. Proceedings of IEEE, 77-2, p. 267-295.

Kaufman, L. and Rousseeuw, P.J. (1990). Finding Groups in Data: An Introduction to Cluster Analysis. Wiley, New York.

Anja Struyf, Mia Hubert and Peter J. Rousseeuw (1996) Clustering in an Object-Oriented Environment. Journal of Statistical Software 1. http://www.jstatsoft.org/v01/i04

Struyf, A., Hubert, M. and Rousseeuw, P.J. (1997). Integrating Robust Clustering Techniques in S-PLUS, Computational Statistics and Data Analysis, 26, 17–37.

Lance, G.N., and W.T. Williams (1966). A General Theory of Classifactory Sorting Strategies, I. Hierarchical Systems. Computer J. 9, 373–380.

Belbin, L., Faith, D.P. and Milligan, G.W. (1992). A Comparison of Two Approaches to Beta-Flexible Clustering. Multivariate Behavioral Research, 27, 417–433.

Gower, J. C. (1971) A general coefficient of similarity and some of its properties, Biometrics 27, 857–874.

Garriga, J., Palmer, J.R., Oltra, A. and Bartumeus, F., 2016. Expectation-maximization binary clustering for behavioural annotation. PLoS One, 11(3), p.e0151984.

Garriga, J., Palmer, J.R.B., Oltra, A. and Bartumeus, F., 2014. EMbC: expectation-maximization binary clustering. arXiv preprint arxiv:1503.04059, 1.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/swift.R
\docType{data}
\name{swift}
\alias{swift}
\title{PAM data collected from Alpine Swift}
\format{
List of 6 variables
\describe{
\item{id}{Logger id}
\item{pressure}{Date in UTC and Atmospheric pressure measurements in hectopascals (every 15 minutes)}
\item{light}{Date in UTC and Light intensity (every 2 minutes)}
\item{acceleration}{Date in UTC, pitch of the accelerometer (i.e. an angle), and derived activity (every 5 minutes)}
\item{temperature}{Date in UTC and temperature on bird in celcius(every 15 minutes)}
\item{magnetic}{Date in UTC and magnetic measurements on 3 axes(every 4 hours)}
}
}
\source{
<www.vogelwarte.ch/en>
}
\usage{
swift
}
\description{
A dataset containing pressure, light, acceleration, temperature and magnetism from an alpine swift (id = 16FY)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/classification2PAM.R
\name{classification2PAM}
\alias{classification2PAM}
\title{Timetable to timeseries}
\usage{
classification2PAM(from, to, classification, addTO, missing = NA)
}
\arguments{
\item{from}{start of event that was classified (generally SOARprep output)}

\item{to}{end of event that was classified (generally SOARprep output)}

\item{classification}{of event (generally classifyPAM()$states output )}

\item{addTO}{data which the classifications are to be added to (e.g. PAM_data$pressure)}

\item{missing}{Missing value replacement. By default NA.}
}
\value{
the classification in addTO dataset
}
\description{
convert a classification timetable into a classification timeseries
}
\examples{
#data(bee_eater)
#PAM_data = bee_eater

#twl = GeoLight::twilightCalc(PAM_data$light$date, PAM_data$light$obs,
#                             LightThreshold = 2, ask = FALSE)


#TOclassify = pamPREP(PAM_data,
#                     method="pressure",
#                     twl = twl,
#                     Pdiff_thld = 2,
#                     light_thld = 2)

#classification = classifyPAM((TOclassify$total_daily_duration *
#                              log(TOclassify$night_P_diff+0.001 )
#                              * TOclassify$total_daily_P_change),
#                             states=3, "hmm")$cluster

#pressure_classification = classification2PAM(from = TOclassify$start,
#                                             to =TOclassify$end,
#                                             classification = classification,
#                                             addTO = PAM_data$pressure,
#                                             missing = NA)

#pressure_classification[pressure_classification == NA] = 0

#plot(PAM_data$pressure$date, PAM_data$pressure$obs,
#     col= viridis::viridis(4)[pressure_classification+1],
#     type="o", pch=16, cex=0.6)

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/quickPLOT.R
\name{quickPLOT}
\alias{quickPLOT}
\title{Quickly plot the pam data to have an idea of it's quality}
\usage{
quickPLOT(
  dta,
  measurements = c("pressure", "light", "acceleration", "temperature", "magnetic"),
  ...
)
}
\arguments{
\item{dta}{path where files are stored}

\item{measurements}{a series of measurements logged by the PAM logger which are to be plotted. Currently supports these file extentions: "pressure","light", "acceleration", "temperature" and "magnetic"}

\item{...}{any additional parameters used by graphics::plot}
}
\value{
a plot of PAM data
}
\description{
Quickly plot the pam data to have an idea of it's quality
}
\examples{
#PAM_data = hoopoe

##plot everything in 2 windows
#par(mar=c(2.5,4,0.5,1))
#quickPLOT(PAM_data)

## only subset some measurements
#quickPLOT(PAM_data, measurements = c("light", "pressure", "acceleration"))

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calculate_water_pressure.R
\name{calculate_water_pressure}
\alias{calculate_water_pressure}
\title{Calculate water pressure from depth}
\usage{
calculate_water_pressure(h, p = 1023.6, g = 9.80665)
}
\arguments{
\item{h}{depth in m *must be positive}

\item{p}{water density for seawater is 1023.6kg/m3 and freshwater 997.0474 kg/m3}

\item{g}{gravitational field strength default 9.80665 m/s2}
}
\value{
pressure in hectopascals
}
\description{
This function calculates water pressure from depth
}
\examples{
pressure = calculate_water_pressure(113.16)

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/compare_classifications.R
\name{compare_classifications}
\alias{compare_classifications}
\title{Combine classification results}
\usage{
compare_classifications(date, classifications)
}
\arguments{
\item{date}{datetime in POSIXCT format. See hoopoe$activity$date for example}

\item{classifications}{a dataframe containing the results from different classifications in each column. The classifications must correspond to the same datetimes.}
}
\value{
a dataframe containing date

the raw classifications

each state with the number of times it was used in a classification

whether or not all classes provided the same state
}
\description{
this function takes multiple classifications of the same data and compares the agreement between each, on a timepoint level.
}
\examples{

# Import data
data(hoopoe)
start = as.POSIXct("2016-07-01","\%Y-\%m-\%d", tz="UTC")
end = as.POSIXct("2017-06-01","\%Y-\%m-\%d", tz="UTC")
PAM_data= create_crop(hoopoe,start,end)

# perform one classification using classifyFLAP
classification = classify_flap(dta = PAM_data$acceleration, period = 12)
# Put the classification in the same resolution as pressure
class1 = create_merged_classification(from = classification$timetable$start,
                            to = classification$timetable$end,
                            # because the timetable only contains migration periods
                            classification = rep_len(1,length(classification$timetable$end)),
                            add_to = PAM_data$pressure)
# Convert to categories
class1 = ifelse(class1 == classification$migration, "Migration", "Other")


# Perform another classification using pressure difference
class2 = c(0,ifelse(abs(diff(PAM_data$pressure$obs))>2, "Migration", "Other"))

# both classes have been converted to the same time intervals as pressure,
date = PAM_data$pressure$date

# Combine the classifications into a dataframe
classifications = data.frame(flap = class1, # flapping classification
                            Pdiff = class2) # pressure difference classification

class_comparison = compare_classifications(date=date,
                                classifications=classifications)

plot(PAM_data$pressure$date,
     PAM_data$pressure$obs,
     type = "l",
     xlab = "Date",
     ylab = "Pressure (hPa)",
     col = "royalblue3")

points(PAM_data$pressure$date,
       PAM_data$pressure$obs,
       cex = class_comparison$Migration / 2,
       col ="orange",
       pch = 16)

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_sensorimage_twilight.R
\name{as_hour}
\alias{as_hour}
\alias{hour_offset}
\title{Hour Manipulation 2}
\usage{
as_hour(tm)

hour_offset(hr, offset = 0)
}
\arguments{
\item{tm}{them timestamp as POSIXct.}

\item{hr}{the decimal hour to be wrap.}

\item{offset}{minimum hour of the interval to wrap into.}
}
\value{
Return a decimal hour.
}
\description{
Utilities for manipulating hours
}
\details{
Given a vector of POSIXct dates, \code{as.hour} extracts the time
of day component of the date and returns it as decimal hours.
Given a vector of decimal hours, \code{hourOffset} recodes the
decimal hour into a new 24 hour interval.
}
\examples{
as_hour(as.POSIXct("2005-11-12 19:58:00"))
hour_offset(1:10,5)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pam3D.R
\name{pam3D}
\alias{pam3D}
\title{3d scatterplot}
\usage{
pam3D(x, y, z, ...)
}
\arguments{
\item{x}{data to plot on x axis}

\item{y}{data to plot on y axis}

\item{z}{data to plot on z axis}

\item{...}{any additional parameters used by rgl::plot3d}
}
\value{
a 3d scatter plot
}
\description{
Creates an interactive 3d scatterplot
}
\examples{
#data("swift")
#start = as.POSIXct("2016-09-01","\%Y-\%m-\%d", tz="UTC")
#end = as.POSIXct("2017-04-15","\%Y-\%m-\%d", tz="UTC")
#swift = cutPAM(swift, start, end)
#PAM_data = swift

#calibration = triMAG(dta = PAM_data$magnetic)

#pam3D(PAM_data$magnetic$mX, PAM_data$magnetic$mY, PAM_data$magnetic$mZ,
#       xlab= "X", ylab= "Y", zlab= "Z")

}
\references{
Adler, D., Nenadic, O. and Zucchini, W., 2003, March. Rgl: A r-library for 3d visualization with opengl. In Proceedings of the 35th Symposium of the Interface: Computing Science and Statistics, Salt Lake City (Vol. 35).
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/changePAM.R
\name{changePAM}
\alias{changePAM}
\title{Changepoint analysis}
\usage{
changePAM(
  dta,
  cpt.method = "meanvar",
  method = "PELT",
  penalty = "Manual",
  pen.value = "100*log(n)",
  ...
)
}
\arguments{
\item{dta}{data to be classified (vector)}

\item{cpt.method}{method for classifying data, currently support "mean","variance", "meanvar", use ?changepoint::cpt.mean ?changepoint::cpt.var or ?changepoint::cpt.meanvar for additional parameters}

\item{method}{changepoint method used in package changepoint. see changepoint::cpt.mean, changepoint::cpt.var and changepoint::cpt.meanvar for details}

\item{penalty}{See changepoint  package for details. Choice of "None", "SIC", "BIC", "MBIC", AIC", "Hannan-Quinn", "Asymptotic", "Manual" and "CROPS" penalties. If Manual is specified, the manual penalty is contained in the pen.value parameter. If Asymptotic is specified, the theoretical type I error is contained in the pen.value parameter. If CROPS is specified, the penalty range is contained in the pen.value parameter; note this is a vector of length 2 which contains the minimum and maximum penalty value. Note CROPS can only be used if the method is "PELT". The predefined penalties listed DO count the changepoint as a parameter, postfix a 0 e.g."SIC0" to NOT count the changepoint as a parameter.}

\item{pen.value}{See changepoint package for details. Tthe theoretical type I error e.g.0.05 when using the Asymptotic penalty. A vector of length 2 (min,max) if using the CROPS penalty. The value of the penalty when using the Manual penalty option - this can be a numeric value or text giving the formula to use. Available variables are, n=length of original data, null=null likelihood, alt=alternative likelihood, tau=proposed changepoint, diffparam=difference in number of alternatve and null parameters}

\item{...}{any additional inputs for changepoint::cpt.mean,  changepoint::cpt.var or changepoint::cpt.meanvar}
}
\value{
Changepoints in the data.
}
\description{
This function implements a changepoint analysis to find changes in variance, mean or both. By default it is set to look for changes in the pressure data to identify migration periods, however it can be altered to any form of changepoint analysis in the package changepoint.
}
\examples{

## Import and crop PAM data
#data("swift")
#start = as.POSIXct("2016-09-01","\%Y-\%m-\%d", tz="UTC")
#end = as.POSIXct("2017-04-15","\%Y-\%m-\%d", tz="UTC")
#swift = cutPAM(swift, start, end)
#PAM_data = swift

# data(bee_eater)
# start = as.POSIXct("2015-08-01","\%Y-\%m-\%d", tz="UTC")
# end = as.POSIXct("2016-06-01","\%Y-\%m-\%d", tz="UTC")
# PAM_data = cutPAM(bee_eater, start, end)

#changepoints  = changePAM(PAM_data$pressure$obs)

# plot the timeseries with the changepoint
#par(mfrow=c(2,1))
#plot(PAM_data$pressure$obs, type="l")
#abline(v=changepoints$changepoints, col="red",lwd=2)

## plot using changepoint package output
#changepoint::plot(changepoints$output, cpt.width=3)

}
\references{
Killick, R. and Eckley, I., 2014. changepoint: An R package for changepoint analysis. Journal of statistical software, 58(3), pp.1-19.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calculate_clockdrift.R
\name{calculate_clockdrift}
\alias{calculate_clockdrift}
\title{Clock Drift Adjustment}
\usage{
calculate_clockdrift(time, start, end)
}
\arguments{
\item{time}{a vector of POSIXct times.}

\item{start}{new start time as POSIXct.}

\item{end}{new end time as POSIXct.}
}
\description{
Adjust time for clock drift
}
\details{
Linearly rescale a sequence of dates to a new start and end time
to correct for clock drift in the tag.
}
\examples{
dates = hoopoe$magnetic$date
drift = 12
adjusted_end = as.POSIXct(dplyr::last(hoopoe$magnetic$date) + 12*60, tz="UTC")
drift_corrected = calculate_clockdrift(time = dates,
                             start= dates[1],
                             end = adjusted_end)

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hoopoe.R
\docType{data}
\name{hoopoe}
\alias{hoopoe}
\title{PAM data collected from Hoopoe}
\format{
List of 6 variables
\describe{
\item{id}{Logger id}
\item{pressure}{Date in UTC and Atmospheric pressure measurements in hectopascals (every 15 minutes)}
\item{light}{Date in UTC and Light intensity (every 5 minutes)}
\item{acceleration}{Date in UTC, pitch of the accelerometer (i.e. an angle), and derived activity (every 5 minutes)}
\item{temperature}{Date in UTC and temperature on bird in celcius(every 15 minutes)}
\item{magnetic}{Date in UTC and magnetic measurements on 3 axes(every 15 minutes)}
}
}
\source{
<www.vogelwarte.ch/en>
}
\usage{
hoopoe
}
\description{
A dataset containint pressure, light, acceleration, temperature and magnetism from a hoopoe (id = 16AJ)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pamSPHERE.R
\name{pamSPHERE}
\alias{pamSPHERE}
\title{Plot 3d sphere}
\usage{
pamSPHERE(
  x,
  y,
  z,
  spherecolor = "white",
  linecolor = "black",
  linewidth = 2,
  ptcol = "black",
  ptsize = 0.02,
  arrows = TRUE,
  cex = 1.5,
  ...
)
}
\arguments{
\item{x}{data on x axis}

\item{y}{data on y axis}

\item{z}{data on z axis}

\item{spherecolor}{color of the sphere}

\item{linecolor}{color of the lines on the sphere}

\item{linewidth}{width of the lines on the sphere}

\item{ptcol}{color of points}

\item{ptsize}{size of the points plotted on the sphere}

\item{arrows}{default TRUE - whether or not to draw the arrows}

\item{cex}{size of font. If cex = 2 then fond is two times bigger.}

\item{...}{additional input for rgl::spheres3d for plotting points on the sphere using}
}
\value{
a 3d sphere plot
}
\description{
This function takes a spherical projection of acceleration or magentic data and plots it on a sphere which can be turned and zoomed into.
}
\examples{
#data("swift")
#start = as.POSIXct("2016-09-01","\%Y-\%m-\%d", tz="UTC")
#end = as.POSIXct("2017-04-15","\%Y-\%m-\%d", tz="UTC")
#swift = cutPAM(swift, start, end)
#PAM_data = swift

## plot an m-phere
#calibration = triMAG(dta = PAM_data$magnetic)
#pamSPHERE(x = calibration$calib_magx,
#          y = calibration$calib_magy,
#          z = calibration$calib_magz,
#          ptcol = "goldenrod",
#          ptsize = 0.03,
#          linecolor ="black",
#          spherecolor="royalblue4",
#          arrows=TRUE,
#          cex=2)

## plot an g-phere
#calibration = triACC(dta = PAM_data$magnetic)
#pamSPHERE(x = calibration$centered_accx,
#          y = calibration$centered_accy,
#          z = calibration$centered_accz,
#          ptcol = "royalblue4",
#          ptsize = 0.03,
#          linecolor ="black",
#          spherecolor="goldenrod",
#          arrows=TRUE)

}
\references{
Adler, D., Nenadic, O. and Zucchini, W., 2003, March. Rgl: A r-library for 3d visualization with opengl. In Proceedings of the 35th Symposium of the Interface: Computing Science and Statistics, Salt Lake City (Vol. 35).
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calculate_depth.R
\name{calculate_depth}
\alias{calculate_depth}
\title{Calculate depth from pressure}
\usage{
calculate_depth(P, p = 1023.6, g = 9.80665)
}
\arguments{
\item{P}{pressure in hectopascals}

\item{p}{water density for seawater is 1023.6kg/m3 and freshwater 997.0474 kg/m3}

\item{g}{gravitational field strength default 9.80665 m/s2}
}
\value{
depth in m
}
\description{
This function calculates depth from pressure
}
\examples{
depth = calculate_depth(11065)

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/create_custom_interpolation.R
\name{create_custom_interpolation}
\alias{create_custom_interpolation}
\title{Reduces the data to a specified temporal resolution}
\usage{
create_custom_interpolation(dta, varint, interp = TRUE, summary = "median")
}
\arguments{
\item{dta}{raw pam data see \verb{data(bee_eater} for example}

\item{varint}{the variable of interest.Sipprots c("pressure","light","acceleration","temperature","magnetic"). For instance if all other variables should be summarised at the same temporal resolution as this specified variable.}

\item{interp}{Default TRUE. Whether or not to replace NAs with interpolated values}

\item{summary}{Can be "sum", "median" or "none". What type of summary variable to give when condensing data -}
}
\value{
reduced/summarised and interpolated dataset
}
\description{
This function takesthe typical PAM_data input, which is a nested list of different sensor data, all formatted at different time resolutions. By specifying a particular sensor, all data are formatted at the same temporal resolution as this data. By default, the median of data with smaller temporal resolution are kept, and data are interpolated.
}
\examples{
data(bee_eater)
PAM_data = bee_eater

start = as.POSIXct("2015-08-01","\%Y-\%m-\%d", tz="UTC")
end = as.POSIXct("2016-06-21","\%Y-\%m-\%d", tz="UTC")
PAM_data = create_crop(PAM_data, start, end)

reduced_dta = create_custom_interpolation(PAM_data , "pressure", interp = FALSE)
head(reduced_dta)

reduced_dta = create_custom_interpolation(PAM_data , "act", interp = FALSE)
head(reduced_dta)

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_sensorimage_legend.R
\name{plot_sensorimage_legend}
\alias{plot_sensorimage_legend}
\title{sensor image legend}
\usage{
plot_sensorimage_legend(
  dta,
  col_palette,
  inset = c(0.5, -0.12),
  cex = 1,
  ncol = 5,
  bg = NA,
  ...
)
}
\arguments{
\item{dta}{data to be plotted (for instance temperature or pressure)}

\item{col_palette}{color palette to use in the graphic}

\item{inset}{where to place the legend default is at the bottom -0.12}

\item{cex}{size of labels}

\item{ncol}{number of columns in legend, default is 5}

\item{bg}{whether or not to have a background, default is NA}

\item{...}{Any additional parameters used by graphics::legend}
}
\value{
an image of the sensor data, for instance with activity it would produce an actogram
}
\description{
This function adds a legend to a sensor image
}
\examples{
#specify the data location
data(hoopoe)
start = as.POSIXct("2016-07-01","\%Y-\%m-\%d", tz="UTC")
end = as.POSIXct("2017-06-01","\%Y-\%m-\%d", tz="UTC")
PAM_data = create_crop(hoopoe,start,end)

# make margins big so that there is enough space for margin
par( mfrow= c(1,1), oma=c(4,2,0,6))

col_palette = c("black",viridis::cividis(90))
par(mar =  c(4,2,4,2))
plot_sensorimage(PAM_data$acceleration$date,
          PAM_data$acceleration$act, main = "Activity",
          col=col_palette, cex=1.2, cex.main = 2)
plot_sensorimage_legend(PAM_data$acceleration$act, col_palette)

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/classify_summary_statistics.R
\name{classify_summary_statistics}
\alias{classify_summary_statistics}
\title{Derive classification data for soaring birds}
\usage{
classify_summary_statistics(
  dta,
  method = "hmm",
  states = 2,
  family = stats::gaussian(),
  ...
)
}
\arguments{
\item{dta}{data to be classified, can be anything}

\item{method}{method for classifying data, currently support "hmm" for hidden markov model and "kmeans" for kmeans clustering}

\item{states}{number of states to classify the data into}

\item{family}{By default \code{gaussian()}. Parameter only used in the hmm method. Must be the same as the parameters used by \code{depmixS4::depmix}}

\item{...}{any additional inputs for depmixs4::depmix, stats::kmeans, cluster::diana, cluster::diana  or EMbC::embc functions, depending on method selected}
}
\value{
the data's classification based on the chosen algorithm
}
\description{
This function takes data, typically a dataframe of events, which are classified into different behavioural states using different algorithms.
}
\examples{

#######################################################
#  data prep
#######################################################

data(bee_eater)
start = as.POSIXct("2015-07-01","\%Y-\%m-\%d", tz="UTC")
end = as.POSIXct("2016-06-01","\%Y-\%m-\%d", tz="UTC")
PAM_data = create_crop(bee_eater, start, end)

twl = GeoLight::twilightCalc(PAM_data$light$date, PAM_data$light$obs,
LightThreshold = 2, ask = FALSE)
availavariable = c("pressure", "light", "acceleration")

to_classify= create_summary_statistics(PAM_data,
                      method= "flap",
                      twl = twl)

to_classify= to_classify[complete.cases(to_classify),]

#######################################################
# k-means example
#######################################################

classification = classify_summary_statistics(to_classify[,c("cum_altitude_change",
                             "night_P_diff" )],
                             states=2, "kmeans")$cluster
pressure_classification = create_merged_classification(from = to_classify$start,
to =to_classify$end,
classification = classification,
add_to = PAM_data$pressure)

plot(PAM_data$pressure$date, PAM_data$pressure$obs,
     type="l")
points(PAM_data$pressure$date, PAM_data$pressure$obs,
       col= pressure_classification+1,
       pch=16)

#######################################################
# HMM example
#######################################################

classification = classify_summary_statistics(to_classify[,c("cum_altitude_change",
                             "night_P_diff" )]
                             #to_classify$night_P_diff +
                             #to_classify$cum_altitude_change
                             ,
                             states=2, "hmm")$cluster
pressure_classification = create_merged_classification(from = to_classify$start,
to =to_classify$end,
classification = classification,
add_to = PAM_data$pressure)

plot(PAM_data$pressure$date, PAM_data$pressure$obs,
     type="l")
points(PAM_data$pressure$date, PAM_data$pressure$obs,
       col= pressure_classification+1,
       pch=16)


#######################################################
# EMBC example
#######################################################

classification = classify_summary_statistics(to_classify[,c("cum_altitude_change",
                             "night_P_diff" )],
                             "embc")

pressure_classification = create_merged_classification(from = to_classify$start,
                                             to =to_classify$end,
                                             classification = classification$cluster,
                                             add_to = PAM_data$pressure)

plot(PAM_data$pressure$date, PAM_data$pressure$obs,
     type="l")
points(PAM_data$pressure$date, PAM_data$pressure$obs,
       col= classification$output@C[pressure_classification+1],
       pch=16)


#######################################################
# agnes example
#######################################################

classification = classify_summary_statistics(to_classify[,c("cum_altitude_change",
                             "night_P_diff" )],
                             states = 2,
                             "agnes")

plot(classification$output, main="agnes", which.plot = 2)

pressure_classification = create_merged_classification(from = to_classify$start,
                                             to =to_classify$end,
                                             classification = classification$cluster,
                                             add_to = PAM_data$pressure)

plot(PAM_data$pressure$date, PAM_data$pressure$obs,
     type="l")
points(PAM_data$pressure$date, PAM_data$pressure$obs,
       col= pressure_classification+1,
       pch=16)

#######################################################
# diana example
#######################################################

classification = classify_summary_statistics(to_classify[,c("cum_altitude_change",
                             "night_P_diff" )],
                             states = 2,
                             "diana")

plot(classification$output, which.plot = 2, main="diana")

pressure_classification = create_merged_classification(from = to_classify$start,
                                             to =to_classify$end,
                                             classification = classification$cluster,
                                             add_to = PAM_data$pressure)

plot(PAM_data$pressure$date, PAM_data$pressure$obs,
     type="l")
points(PAM_data$pressure$date, PAM_data$pressure$obs,
       col= pressure_classification+1,
       pch=16)



}
\references{
Forgy, E. W. (1965). Cluster analysis of multivariate data: efficiency vs interpretability of classifications. Biometrics, 21, 768–769.

Hartigan, J. A. and Wong, M. A. (1979). Algorithm AS 136: A K-means clustering algorithm. Applied Statistics, 28, 100–108. doi: 10.2307/2346830.

Lloyd, S. P. (1957, 1982). Least squares quantization in PCM. Technical Note, Bell Laboratories. Published in 1982 in IEEE Transactions on Information Theory, 28, 128–137.

MacQueen, J. (1967). Some methods for classification and analysis of multivariate observations. In Proceedings of the Fifth Berkeley Symposium on Mathematical Statistics and Probability, eds L. M. Le Cam & J. Neyman, 1, pp. 281–297. Berkeley, CA: University of California Press.

Ingmar Visser and Maarten Speekenbrink (2010). depmixS4: An R Package for Hidden Markov Models. Journal of Statistical Software, 36(7), p. 1-21.

Lawrence R. Rabiner (1989). A tutorial on hidden Markov models and selected applications in speech recognition. Proceedings of IEEE, 77-2, p. 267-295.

Kaufman, L. and Rousseeuw, P.J. (1990). Finding Groups in Data: An Introduction to Cluster Analysis. Wiley, New York.

Anja Struyf, Mia Hubert and Peter J. Rousseeuw (1996) Clustering in an Object-Oriented Environment. Journal of Statistical Software 1. http://www.jstatsoft.org/v01/i04

Struyf, A., Hubert, M. and Rousseeuw, P.J. (1997). Integrating Robust Clustering Techniques in S-PLUS, Computational Statistics and Data Analysis, 26, 17–37.

Lance, G.N., and W.T. Williams (1966). A General Theory of Classifactory Sorting Strategies, I. Hierarchical Systems. Computer J. 9, 373–380.

Belbin, L., Faith, D.P. and Milligan, G.W. (1992). A Comparison of Two Approaches to Beta-Flexible Clustering. Multivariate Behavioral Research, 27, 417–433.

Gower, J. C. (1971) A general coefficient of similarity and some of its properties, Biometrics 27, 857–874.

Garriga, J., Palmer, J.R., Oltra, A. and Bartumeus, F., 2016. Expectation-maximization binary clustering for behavioural annotation. PLoS One, 11(3), p.e0151984.

Garriga, J., Palmer, J.R.B., Oltra, A. and Bartumeus, F., 2014. EMbC: expectation-maximization binary clustering. arXiv preprint arxiv:1503.04059, 1.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cutPAM.R
\name{cutPAM}
\alias{cutPAM}
\title{Crop all sensor data to the same time}
\usage{
cutPAM(dta, start, end)
}
\arguments{
\item{dta}{path where files are stored}

\item{start}{posicxt object for date that PAM data should start}

\item{end}{posicxt object for date that PAM data should end}
}
\value{
shortened PAM data
}
\description{
Get rid of excess data. e.g. when a logger is kept in a rucksack or a lab before being downloaded.
}
\examples{
#data(hoopoe)
#PAM_data=hoopoe
#str(PAM_data)

#start = as.POSIXct("2016-07-01","\%Y-\%m-\%d", tz="UTC")
#end = as.POSIXct("2017-06-01","\%Y-\%m-\%d", tz="UTC")

#newPAM = cutPAM(PAM_data,start,end)

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/create_merge.R
\name{create_merge}
\alias{create_merge}
\title{Merges all sensors into one}
\usage{
create_merge(dta, all = TRUE, interp = FALSE)
}
\arguments{
\item{dta}{raw pam data see \verb{data(bee_eater} for example}

\item{all}{logical. Default TRUE. Whether or not to keep NAs (i.e. all the datasets)}

\item{interp}{logical. Default FALSE. whether or not to interpolate if there are NAs. if all = FALSE then interp is not used.}
}
\value{
merged and interpolated dataset
}
\description{
This function takesthe typical PAM_data input, which is a nested list of different sensor data, all formatted at different time resolutions, and merges them all into one big table. By default all times are kept, and not interpolated.
}
\examples{
data(bee_eater)
PAM_data = bee_eater
merged_dta = create_merge(PAM_data, interp = TRUE)
head(merged_dta)

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bee_eater.R
\docType{data}
\name{bee_eater}
\alias{bee_eater}
\title{PAM data collected from European Bee-Eater}
\format{
List of 6 variables
\describe{
\item{id}{Logger id}
\item{pressure}{Date in UTC and Atmospheric pressure measurements in hectopascals (every 15 minutes)}
\item{light}{Date in UTC and Light intensity (every 5 minutes)}
\item{acceleration}{Date in UTC, pitch of the accelerometer (i.e. an angle), and derived activity (every 5 minutes)}
\item{temperature}{Date in UTC and temperature on bird in celcius(every 15 minutes)}
\item{magnetic}{Date in UTC and magnetic measurements on 3 axes(every 15 minutes)}
}
}
\source{
<www.vogelwarte.ch/en>
}
\usage{
bee_eater
}
\description{
A dataset containint pressure, light, acceleration, temperature and magnetism from a european bee-eater (id = 14SA)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_sensorimage_twilight.R
\name{plot_sensorimage_twilight}
\alias{plot_sensorimage_twilight}
\title{Add sunrise and unset to a sensor image}
\usage{
plot_sensorimage_twilight(date, offset, ...)
}
\arguments{
\item{date}{Date data in POSIXct format, most commonly \code{PAM_data$acceleration$date}}

\item{offset}{This parameter determines where the center of the graph is. When \code{offset = 0}, then midday is at the center of the graph. when \code{offset=12} midnight`}

\item{...}{Any additional parameters taken by graphics::points which the user may want to use to modify the graphic}
}
\value{
a plot
}
\description{
Add sunrise and unset to a sensor image
}
\examples{

data(hoopoe)

# Calculate twilights
twilights <- GeoLight::twilightCalc(hoopoe$light$date, hoopoe$light$obs,
LightThreshold = 2, ask = FALSE)

par(mar=c(4,4,4,2),mfrow=c(1,2), oma= c(0,0,0,8))

# Plot daytime in middle
plot_sensorimage(hoopoe$acceleration$date,
          hoopoe$acceleration$act,
          main = "Day in middle",ploty=FALSE,
          offset=0,
          col=c("black",viridis::cividis(90)),
          cex=1.2, cex.main = 2)
#Add sunrises and sunsets
plot_sensorimage_twilight(twilights$tFirst,
       offset=0,
       col= ifelse(twilights$type == 1,
                   "goldenrod","cornflowerblue"),
       pch=16, cex=0.5)

# Plot nightime in the middle
offset=12
plot_sensorimage(hoopoe$acceleration$date,
          hoopoe$acceleration$act,
          main = "Night in middle",ploty=TRUE,
          offset=offset,
          col=c("black",viridis::cividis(90)),
          cex=1.2, cex.main = 2)

plot_sensorimage_twilight(twilights$tFirst,
       offset=offset,
       col= ifelse(twilights$type == 1,
                   "goldenrod","cornflowerblue"),
       pch=16, cex=0.5)



}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dygraphClassified.R
\name{dygraphClassified}
\alias{dygraphClassified}
\title{Plot PAM data with dygraphs}
\usage{
dygraphClassified(
  dta,
  from = dta$light$date[1],
  to = dta$light$date[length(dta$light$date)],
  toPLOT = names(dta),
  timetable = timetable
)
}
\arguments{
\item{dta}{PAM data to be plotted}

\item{from}{date that plotting starts}

\item{to}{date that plotting ends}

\item{toPLOT}{names of the variables to plot. For now this includes \code{light}, \code{pressure}, \code{acceleration} and \code{temperature}}

\item{timetable}{classification of start/stop}
}
\value{
a plot of all the measurements
}
\description{
This opens a java application which allows the user to zoom in and out. In Rstudio it will open in the viewer pane and in base R in an html readers. Note that this can be a bit slow
}
\examples{
##load dummy data
#data(hoopoe)
#PAM_data=hoopoe

## This bit is for Rstudio users to prevent html from opening in Viewer pane and crashing
## It opens in web browser instead
#backup_options <- options()
#options(viewer=NULL)

## Classify bir's behaviour
#behaviour = classifyFLAP(dta = PAM_data$acceleration, period = 5)

## Plot interactive graphics
#dygraphClassified(dta = PAM_data, timetable = behaviour$timetable)

## restore Rstudio settings from before plot
#options(backup_options)

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_interactive_timeseries_classification.R
\name{plot_interactive_timeseries_classification}
\alias{plot_interactive_timeseries_classification}
\title{Plot PAM data as an interactive timeseries}
\usage{
plot_interactive_timeseries_classification(
  to_classify,
  classification_datetime,
  classification,
  to_plot = c("light", "pressure", "temperature", "acceleration", "magnetic")
)
}
\arguments{
\item{to_classify}{rolling window data to plot}

\item{classification_datetime}{dates associated with classification}

\item{classification}{output of the classification}

\item{to_plot}{names of the variables to plot. For now this includes \code{light}, \code{pressure}, \code{acceleration} and \code{temperature}}
}
\value{
a plot of all the measurements
}
\description{
This opens a java dygraph application which allows the user to zoom in and out. In Rstudio it will open in the viewer pane and in base R in an html readers. Note that this can be a bit slow
}
\examples{
#load dummy data
data(hoopoe)
PAM_data=hoopoe

# This bit is for Rstudio users to prevent html from opening in Viewer pane and crashing
# It opens in web browser instead
backup_options <- options()
options(viewer=NULL)

# Plot interactive graphics
plot_interactive_timeseries(dta = PAM_data)

# restore Rstudio settings from before plot
options(backup_options)


}
\references{
Vanderkam, D., Allaire, J., Owen, J., Gromer, D., Shevtsov, P. and Thieurmel, B., dygraphs: Interface to Dygraphs Interactive Time Series Charting Library, 2015. URL http://CRAN. R-project. org/package= dygraphs. R package version 0.4, 5, p.7.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/triMAG.R
\name{triMAG}
\alias{triMAG}
\title{Calibrate magnetic data}
\usage{
triMAG(dta)
}
\arguments{
\item{dta}{magentic data from PAM logger}
}
\value{
roll, pitch and yaw from acceleration data, as well as calibrated magnetic data for x, y and z axes
}
\description{
This function calibrates tri-axial magnetic data and then also calculates yaw, pitch and roll
}
\examples{
#data(swift)
#PAM_data = swift

#calibration = triMAG(dta = PAM_data$magnetic)

}
\references{
Bidder, O.R., Walker, J.S., Jones, M.W., Holton, M.D., Urge, P., Scantlebury, D.M., Marks, N.J., Magowan, E.A., Maguire, I.E. and Wilson, R.P., 2015. Step by step: reconstruction of terrestrial animal movement paths by dead-reckoning. Movement ecology, 3(1), p.23.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/classifyFLAP.R
\name{classifyFLAP}
\alias{classifyFLAP}
\title{Classify flapping flight}
\usage{
classifyFLAP(dta, period = 12, toPLOT = TRUE, method = "kmeans", tz = "UTC")
}
\arguments{
\item{dta}{data stored as a list see str(data(PAM_data)) for example format}

\item{period}{number of timepoints after which behaviour is considered migratory e.g. for hoopoes, 3x5min = 15 minutes of intense activity is considered flight}

\item{toPLOT}{can be true or false. If true then threshold is plotted according to plotTHLD()}

\item{method}{for the time being only supports "kmeans", but will later also include maybe}

\item{tz}{timezone, default is "UTC"}
}
\value{
timetable: a timetable for when the species was migrating or not,

classification: a classification timeseries where datetime corresponds to activity, and

no_movement: the value in classification which corresponds to no movement

low_movement: the value in classification which corresponds to low activity

high_movement: the value in classification which corresponds to high activity

migration: the value in classification which corresponds to migratory flapping flight

threshold: the threshold between low and high activity
}
\description{
This function uses activity data to classify migratory flapping flight.
}
\examples{
##specify the data location
#data(hoopoe)
## make sure the cropping period is in the correct date format
#start = as.POSIXct("2016-07-01","\%Y-\%m-\%d", tz="UTC")
#end = as.POSIXct("2017-06-01","\%Y-\%m-\%d", tz="UTC")

## Crop the data
#PAM_data= cutPAM(hoopoe,start,end)
#str(PAM_data)

#behaviour = classifyFLAP(dta = PAM_data$acceleration, period = 12)

#col=c("brown","cyan4","black","gold")
#plot(PAM_data$acceleration$date[6000:8000],PAM_data$acceleration$act[6000:8000],
#     col=col[behaviour$classification][6000:8000],
#     type="o", pch=20,
#     xlab="Date",
#     ylab="Activity")

#behaviour$timetable

}
\references{
Bäckman, J., Andersson, A., Alerstam, T., Pedersen, L., Sjöberg, S., Thorup, K. and Tøttrup, A.P., 2017. Activity and migratory flights of individual free‐flying songbirds throughout the annual cycle: method and first case study. Journal of avian biology, 48(2), pp.309-319.

Liechti, F., Bauer, S., Dhanjal-Adams, K.L., Emmenegger, T., Zehtindjiev, P. and Hahn, S., 2018. Miniaturized multi-sensor loggers provide new insight into year-round flight behaviour of small trans-Sahara avian migrants. Movement ecology, 6(1), p.19.

Bruderer, B., Peter, D., Boldt, A. and Liechti, F., 2010. Wing‐beat characteristics of birds recorded with tracking radar and cine camera. Ibis, 152(2), pp.272-291.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/addTWL.R
\name{addTWL}
\alias{addTWL}
\title{Add sunrise sunset to Actogram}
\usage{
addTWL(date, offset, ...)
}
\arguments{
\item{date}{Date data in POSIXct format, most commonly \code{PAM_data$acceleration$date}}

\item{offset}{This parameter determines where the center of the graph is. When \code{offset = 0}, then midday is at the center of the graph. when \code{offset=12} midnight`}

\item{...}{Any additional parameters taken by graphics::points which the user may want to use to modify the graphic}
}
\value{
a plot
}
\description{
Add sunrise sunset to Actogram
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_timeseries.R
\name{plot_timeseries}
\alias{plot_timeseries}
\title{Quickly plot the pam data to have an idea of it's quality}
\usage{
plot_timeseries(
  dta,
  measurements = c("pressure", "light", "acceleration", "temperature", "magnetic"),
  ...
)
}
\arguments{
\item{dta}{path where files are stored}

\item{measurements}{a series of measurements logged by the PAM logger which are to be plotted. Currently supports these file extentions: "pressure","light", "acceleration", "temperature" and "magnetic"}

\item{...}{any additional parameters used by graphics::plot}
}
\value{
a plot of PAM data
}
\description{
Quickly plot the pam data to have an idea of it's quality
}
\examples{
PAM_data = hoopoe

#plot everything in 2 windows
par(mar=c(2.5,4,0.5,1))
plot_timeseries(PAM_data)

# only subset some measurements
plot_timeseries(PAM_data, measurements = c("light",
      "pressure", "acceleration"))

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calculate_roll.R
\name{calculate_roll}
\alias{calculate_roll}
\alias{calculate_pitch}
\alias{calculate_yaw}
\title{triaxial calculations}
\usage{
calculate_roll(dta)

calculate_pitch(dta)

calculate_yaw(dta)
}
\arguments{
\item{dta}{magentic data from PAM logger see hoopoe$magnetic for an example}
}
\value{
roll, pitch and yaw from acceleration data
}
\description{
Utilities for calculating roll, pitch and yaw
}
\details{
Calculate roll, pitch and yaw
}
\examples{
calculate_roll(dta = swift$magnetic)
calculate_pitch(dta = swift$magnetic)
calculate_yaw(dta = swift$magnetic)

}
\references{
Bidder, O.R., Walker, J.S., Jones, M.W., Holton, M.D., Urge, P., Scantlebury, D.M., Marks, N.J., Magowan, E.A., Maguire, I.E. and Wilson, R.P., 2015. Step by step: reconstruction of terrestrial animal movement paths by dead-reckoning. Movement ecology, 3(1), p.23.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotTHLD.R
\name{plotTHLD}
\alias{plotTHLD}
\title{Plot threshold}
\usage{
plotTHLD(dta, classification, threshold, type, new_window = FALSE, ...)
}
\arguments{
\item{dta}{Raw acceleration or pressure data used to make the classification}

\item{classification}{Is the result of the classification. is in the format of a vector of numbers used for low or high activity}

\item{threshold}{The threshold between different classes}

\item{type}{The type of classification used i.e "flapping" or soar-gliding}

\item{new_window}{whether you want to plot it in a new window or not}

\item{...}{any additional parameters used by graphics::hist}
}
\value{
a graphic with the output from the classifications
}
\description{
Plot threshold
}
\examples{
##specify the data location
#data(hoopoe)
#PAM_data=hoopoe

#PAM_data$acceleration = PAM_data$acceleration[((PAM_data$acceleration$date >= "2016-07-30")
#& (PAM_data$acceleration$date <= "2017-06-01")),]

#behaviour = classifyFLAP(dta = PAM_data$acceleration,
#                         period = 3,
#                         toPLOT = FALSE)

#plotTHLD(dta = PAM_data$acceleration$act,
#         type = "flapping",
#         classification = behaviour$classification,
#         threshold = behaviour$threshold)

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mergePAM.R
\name{mergePAM}
\alias{mergePAM}
\title{Merges all sensors into one}
\usage{
mergePAM(dta, all = TRUE, interp = FALSE)
}
\arguments{
\item{dta}{raw pam data see \verb{data(bee_eater} for example}

\item{all}{logical. Default TRUE. Whether or not to keep NAs (i.e. all the datasets)}

\item{interp}{logical. Default FALSE. whether or not to interpolate if there are NAs. if all = FALSE then interp is not used.}
}
\value{
merged and interpolated dataset
}
\description{
This function takesthe typical PAM_data input, which is a nested list of different sensor data, all formatted at different time resolutions, and merges them all into one big table. By default all times are kept, and not interpolated.
}
\examples{
#data(bee_eater)
#PAM_data = bee_eater
#merged_dta = mergePAM(PAM_data, interp = TRUE)
#head(merged_dta)

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sensorIMG.R
\name{sensorIMG}
\alias{sensorIMG}
\title{sensor image}
\usage{
sensorIMG(
  date,
  sensor_data,
  tz = "UTC",
  plotx = TRUE,
  ploty = TRUE,
  labelx = TRUE,
  labely = TRUE,
  offset = 0,
  dt = NA,
  xlab = "Hour",
  ylab = "Date",
  cex = 2,
  col = c("black", viridis::magma(90)),
  ...
)
}
\arguments{
\item{date}{Date data in POSIXct format, most commonly \code{PAM_data$acceleration$date}}

\item{sensor_data}{sensor data, for example look at \code{PAM_data$acceleration$act}}

\item{tz}{Time zone for POSIXct, default set to "UTC"}

\item{plotx}{wherether or not to plot the x axis ticks + labels (for instance when compiling multifigures)}

\item{ploty}{wherether or not to plot the y axis ticks + labels (for instance when compiling multifigures)}

\item{labelx}{wherether or not to write the name of the x axis (for instance when compiling multifigures)}

\item{labely}{wherether or not to write the name of the y axis (for instance when compiling multifigures)}

\item{offset}{This parameter determines where the center of the graph is. When \code{offset = 0}, then midday is at the center of the graph. when \code{offset=12} midnight`}

\item{dt}{the time interval to which the data are resampled (secs). Default is \code{NA}}

\item{xlab}{label for x-axis (as a character string)}

\item{ylab}{label for y axis (as a character string)}

\item{cex}{size of labels}

\item{col}{Colour scheme of plot. Default \code{col = c("black",viridis::magma(90))}}

\item{...}{Any additional parameters used by graphics::image}
}
\value{
an image of the sensor data, for instance with activity it would produce an actogram
}
\description{
This function plots sensor data as an image. An actogram for instance is a type of sensor image.
}
\examples{
##specify the data location
#data(hoopoe)
#start = as.POSIXct("2016-07-01","\%Y-\%m-\%d", tz="UTC")
#end = as.POSIXct("2017-06-01","\%Y-\%m-\%d", tz="UTC")
#PAM_data = cutPAM(hoopoe,start,end)

## Create plots with 3 together (mfrow)
#par( mfrow= c(1,3), oma=c(0,2,0,6))

#par(mar =  c(4,2,4,2))
#sensorIMG(PAM_data$acceleration$date, ploty=FALSE,
#          PAM_data$acceleration$act, main = "Activity",
#          col=c("black",viridis::cividis(90)), cex=1.2, cex.main = 2)

#par(mar =  c(4,2,4,2))
#sensorIMG(PAM_data$pressure$date, plotx=TRUE, ploty=FALSE, labely=FALSE,
#          PAM_data$pressure$obs,  main="Pressure",
#          col=c("black",viridis::cividis(90)), cex=1.2, cex.main = 2)

#par(mar =  c(4,2,4,2))
#sensorIMG(PAM_data$temperature$date, labely=FALSE,
#          PAM_data$temperature$obs,  main="Temperature",
#          col=c("black",viridis::cividis(90)), cex=1.2, cex.main = 2)

######################################################
# Look at a classification output
######################################################

## Classification
#classification  =  classifyFLAP(dta = PAM_data$acceleration, period = 10, toPLOT=FALSE)

#par( mfrow= c(1,3), oma=c(0,2,0,6),mar =  c(4,2,4,2))

#sensorIMG(PAM_data$pressure$date, c(0,abs(diff(PAM_data$pressure$obs))),
#          main="Pressure  difference",
#          ploty=FALSE,
#          col=c("black",viridis::cividis(90)), cex=1.2, cex.main = 2)

#sensorIMG(PAM_data$acceleration$date, PAM_data$acceleration$act,  main="Activity",
#          ploty=FALSE, labely=FALSE,
#          col=c(viridis::cividis(90)), cex=1.2, cex.main = 2)

#sensorIMG(PAM_data$acceleration$date,
#          ifelse(classification$classification == classification$migration, 1,2),
#          main="Migration Classification",
#          labely=FALSE,
#          col = c("orange","black"),
#          cex=1.2, cex.main = 2)


#twilights <- GeoLight::twilightCalc(PAM_data$light$date,
#                                    PAM_data$light$obs,
#                                    LightThreshold = 2,
#                                    ask = FALSE)

#addTWL(twilights$tFirst, offset=0,
#       col= ifelse(twilights$type == 1,
#                   "goldenrod","cornflowerblue"),
#       pch=16, cex=0.5)


}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/classify_pressurechange.R
\name{classify_pressurechange}
\alias{classify_pressurechange}
\title{Classify pressure change}
\usage{
classify_pressurechange(dta, thld = 2, duration = 1, tz = "UTC")
}
\arguments{
\item{dta}{data stored as a list see str(hoopoe) for example format}

\item{thld}{threshold that is used di discinguish between background pressure fluctuations and flights or dives, by default it is 2hPa}

\item{duration}{duration in hours of the pressure displacement used to quantify a long duration dive or flight, default is 1h}

\item{tz}{timezone, default is "UTC"}
}
\value{
timetable: a timetable for when the species was migrating or not,

classification: a classification timeseries where datetime corresponds to activity, and

no_pressurechange: the value in classification which corresponds to background pressure changes from weather

short_pressurechange: the value in classification which corresponds to a short change in pressure

long_pressurechange: the value in classification which corresponds to a long change in pressure

threshold: the threshold between low and high activity
}
\description{
This function uses a change in pressure greater than a certain threshold to classify diving or flying
}
\examples{
#specify the data location
data(hoopoe)
# make sure the cropping period is in the correct date format
start = as.POSIXct("2016-07-01","\%Y-\%m-\%d", tz="UTC")
end = as.POSIXct("2017-06-01","\%Y-\%m-\%d", tz="UTC")

# Crop the data
PAM_data= create_crop(hoopoe,start,end)

behaviour = classify_pressurechange(dta = PAM_data$pressure)

col=c("black","cyan4","gold")
plot(PAM_data$pressure$date[2000:2800],PAM_data$pressure$obs[2000:2800],
     col=col[behaviour$classification+1][2000:2800],
     type="o", pch=20,
     xlab="Date",
     ylab="Pressure")

behaviour$timetable

}
\references{
Bäckman, J., Andersson, A., Alerstam, T., Pedersen, L., Sjöberg, S., Thorup, K. and Tøttrup, A.P., 2017. Activity and migratory flights of individual free‐flying songbirds throughout the annual cycle: method and first case study. Journal of avian biology, 48(2), pp.309-319.

Liechti, F., Bauer, S., Dhanjal-Adams, K.L., Emmenegger, T., Zehtindjiev, P. and Hahn, S., 2018. Miniaturized multi-sensor loggers provide new insight into year-round flight behaviour of small trans-Sahara avian migrants. Movement ecology, 6(1), p.19.

Bruderer, B., Peter, D., Boldt, A. and Liechti, F., 2010. Wing‐beat characteristics of birds recorded with tracking radar and cine camera. Ibis, 152(2), pp.272-291.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calculate_triaxial_magnetic.R
\name{calculate_triaxial_magnetic}
\alias{calculate_triaxial_magnetic}
\title{Calibrate magnetic data}
\usage{
calculate_triaxial_magnetic(dta)
}
\arguments{
\item{dta}{magentic data from PAM logger}
}
\value{
roll, pitch and yaw from acceleration data, as well as calibrated magnetic data for x, y and z axes
}
\description{
This function calibrates tri-axial magnetic data and then also calculates yaw, pitch and roll
}
\examples{
data(swift)
PAM_data = swift

calibration = calculate_triaxial_magnetic(dta = PAM_data$magnetic)

}
\references{
Bidder, O.R., Walker, J.S., Jones, M.W., Holton, M.D., Urge, P., Scantlebury, D.M., Marks, N.J., Magowan, E.A., Maguire, I.E. and Wilson, R.P., 2015. Step by step: reconstruction of terrestrial animal movement paths by dead-reckoning. Movement ecology, 3(1), p.23.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dygraphPAM.R
\name{dygraphPAM}
\alias{dygraphPAM}
\title{Plot PAM data with dygraphs}
\usage{
dygraphPAM(
  dta,
  from = dta$light$date[1],
  to = dta$light$date[length(dta$light$date)],
  toPLOT = names(dta)
)
}
\arguments{
\item{dta}{PAM data to be plotted}

\item{from}{date that plotting starts}

\item{to}{date that plotting ends}

\item{toPLOT}{names of the variables to plot. For now this includes \code{light}, \code{pressure}, \code{acceleration} and \code{temperature}}
}
\value{
a plot of all the measurements
}
\description{
This opens a java dygraph application which allows the user to zoom in and out. In Rstudio it will open in the viewer pane and in base R in an html readers. Note that this can be a bit slow
}
\examples{
##load dummy data
#data(hoopoe)
#PAM_data=hoopoe

## This bit is for Rstudio users to prevent html from opening in Viewer pane and crashing
## It opens in web browser instead
#backup_options <- options()
#options(viewer=NULL)

## Plot interactive graphics
#dygraphPAM(dta = PAM_data)

## restore Rstudio settings from before plot
#options(backup_options)


}
\references{
Vanderkam, D., Allaire, J., Owen, J., Gromer, D., Shevtsov, P. and Thieurmel, B., dygraphs: Interface to Dygraphs Interactive Time Series Charting Library, 2015. URL http://CRAN. R-project. org/package= dygraphs. R package version 0.4, 5, p.7.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PAMLr.R
\docType{package}
\name{PAMLr}
\alias{PAMLr}
\title{PAMLr}
\description{
This package manipulates data from SOI-GDL3pam loggers (developped by the Swiss Ornithological Institute). These measure Pressure, Activity, Magnetism and Light.
}
\author{
Kiran Dhanjal-Adams \email{kiran.dhanjal.adams@gmail.com}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_interactive_3d.R
\name{plot_interactive_3d}
\alias{plot_interactive_3d}
\title{3d scatterplot}
\usage{
plot_interactive_3d(x, y, z, ...)
}
\arguments{
\item{x}{data to plot on x axis}

\item{y}{data to plot on y axis}

\item{z}{data to plot on z axis}

\item{...}{any additional parameters used by rgl::plot3d}
}
\value{
a 3d scatter plot
}
\description{
Creates an interactive 3d scatterplot
}
\examples{
data("swift")
start = as.POSIXct("2016-09-01","\%Y-\%m-\%d", tz="UTC")
end = as.POSIXct("2017-04-15","\%Y-\%m-\%d", tz="UTC")
swift = create_crop(swift, start, end)
PAM_data = swift

calibration = calculate_triaxial_magnetic(dta = PAM_data$magnetic)

plot_interactive_3d(PAM_data$magnetic$mX, PAM_data$magnetic$mY, PAM_data$magnetic$mZ,
       xlab= "X", ylab= "Y", zlab= "Z")

}
\references{
Adler, D., Nenadic, O. and Zucchini, W., 2003, March. Rgl: A r-library for 3d visualization with opengl. In Proceedings of the 35th Symposium of the Interface: Computing Science and Statistics, Salt Lake City (Vol. 35).
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/classify_flap.R
\name{classify_flap}
\alias{classify_flap}
\title{Classify flapping flight}
\usage{
classify_flap(dta, period = 12, to_plot = TRUE, method = "kmeans", tz = "UTC")
}
\arguments{
\item{dta}{data stored as a list see str(data(PAM_data)) for example format}

\item{period}{number of timepoints after which behaviour is considered migratory e.g. for hoopoes, 3x5min = 15 minutes of intense activity is considered flight}

\item{to_plot}{can be true or false. If true then threshold is plotted according to plot_histogram()}

\item{method}{for the time being only supports "kmeans", but will later also include maybe}

\item{tz}{timezone, default is "UTC"}
}
\value{
timetable: a timetable for when the species was migrating or not,

classification: a classification timeseries where datetime corresponds to activity, and

no_movement: the value in classification which corresponds to no movement

low_movement: the value in classification which corresponds to low activity

high_movement: the value in classification which corresponds to high activity

migration: the value in classification which corresponds to migratory flapping flight

threshold: the threshold between low and high activity
}
\description{
This function uses activity data to classify migratory flapping flight.
}
\examples{
#specify the data location
data(hoopoe)
# make sure the cropping period is in the correct date format
start = as.POSIXct("2016-07-01","\%Y-\%m-\%d", tz="UTC")
end = as.POSIXct("2017-06-01","\%Y-\%m-\%d", tz="UTC")

# Crop the data
PAM_data= create_crop(hoopoe,start,end)
str(PAM_data)

behaviour = classify_flap(dta = PAM_data$acceleration, period = 12)

col=c("brown","cyan4","black","gold")
plot(PAM_data$acceleration$date[6000:8000],PAM_data$acceleration$act[6000:8000],
     col=col[behaviour$classification][6000:8000],
     type="o", pch=20,
     xlab="Date",
     ylab="Activity")

behaviour$timetable

}
\references{
Bäckman, J., Andersson, A., Alerstam, T., Pedersen, L., Sjöberg, S., Thorup, K. and Tøttrup, A.P., 2017. Activity and migratory flights of individual free‐flying songbirds throughout the annual cycle: method and first case study. Journal of avian biology, 48(2), pp.309-319.

Liechti, F., Bauer, S., Dhanjal-Adams, K.L., Emmenegger, T., Zehtindjiev, P. and Hahn, S., 2018. Miniaturized multi-sensor loggers provide new insight into year-round flight behaviour of small trans-Sahara avian migrants. Movement ecology, 6(1), p.19.

Bruderer, B., Peter, D., Boldt, A. and Liechti, F., 2010. Wing‐beat characteristics of birds recorded with tracking radar and cine camera. Ibis, 152(2), pp.272-291.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_histogram.R
\name{plot_histogram}
\alias{plot_histogram}
\title{Plot threshold}
\usage{
plot_histogram(dta, classification, threshold, type, new_window = FALSE, ...)
}
\arguments{
\item{dta}{Raw acceleration or pressure data used to make the classification}

\item{classification}{Is the result of the classification. is in the format of a vector of numbers used for low or high activity}

\item{threshold}{The threshold between different classes}

\item{type}{The type of classification used i.e "flapping" or soar-gliding}

\item{new_window}{whether you want to plot it in a new window or not}

\item{...}{any additional parameters used by graphics::hist}
}
\value{
a graphic with the output from the classifications
}
\description{
Plot threshold
}
\examples{
#specify the data location
data(hoopoe)
PAM_data=hoopoe

PAM_data$acceleration = PAM_data$acceleration[((PAM_data$acceleration$date >= "2016-07-30")
& (PAM_data$acceleration$date <= "2017-06-01")),]

behaviour = classify_flap(dta = PAM_data$acceleration,
                         period = 3,
                         to_plot = FALSE)

plot_histogram(dta = PAM_data$acceleration$act,
         type = "flapping",
         classification = behaviour$classification,
         threshold = behaviour$threshold)

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/clockDRIFT.R
\name{clockDRIFT}
\alias{clockDRIFT}
\title{Clock Drift Adjustment}
\usage{
clockDRIFT(time, start, end)
}
\arguments{
\item{time}{a vector of POSIXct times.}

\item{start}{new start time as POSIXct.}

\item{end}{new end time as POSIXct.}
}
\description{
Adjust time for clock drift
}
\details{
Linearly rescale a sequence of dates to a new start and end time
to correct for clock drift in the tag.
}
\examples{
#dates = hoopoe$magnetic$date
#drift = 12
#adjusted_end = as.POSIXct(dplyr::last(hoopoe$magnetic$date) + 12*60, tz="UTC")
#drift_corrected = clockDRIFT(time = dates,
#                             start= dates[1],
#                             end = adjusted_end)

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/create_crop.R
\name{create_crop}
\alias{create_crop}
\title{Crop all sensor data to the same timeframe}
\usage{
create_crop(dta, start, end)
}
\arguments{
\item{dta}{path where files are stored}

\item{start}{posicxt object for date that PAM data should start}

\item{end}{posicxt object for date that PAM data should end}
}
\value{
shortened PAM data
}
\description{
Get rid of excess data. e.g. when a logger is kept in a rucksack or a lab before being downloaded.
}
\examples{
data(hoopoe)
PAM_data=hoopoe
str(PAM_data)

start = as.POSIXct("2016-07-01","\%Y-\%m-\%d", tz="UTC")
end = as.POSIXct("2017-06-01","\%Y-\%m-\%d", tz="UTC")

newPAM = create_crop(PAM_data,start,end)

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/compare_confusion_matrix.R
\name{compare_confusion_matrix}
\alias{compare_confusion_matrix}
\title{create a confusion matrix}
\usage{
compare_confusion_matrix(reference, classified)
}
\arguments{
\item{reference}{Reference dataset wich classification is to be compared to}

\item{classified}{Classification output to compare with a reference dataset}
}
\value{
a confusion matrix
}
\description{
The function populates a confusion matrix using predicted and reference points and estimate the errors in commission, omissioin for each class, in addition to Producer, User and overall accuracy, and the Kappa Coefficient. Errors in Commission provide a measure of false negatives i.e. the number of points that were predicted to be part of a class that they were not (probability something was incorrectly classified FN/(TP+FN)). Errors in Omission provide a measure of false positives that were predicted to be in a different class from their actual class (probability that something was missed FP/(FP +TP). User Accuracy or Recall represents the probability that a class was correctly classified TP/(TP + FN). Producer Accuracy or Precision provides a measure of how likely something was missed by the classification (probability that something was not missed TP/(TP + FP)). The Overall Accuracy represents the probability that all classes were correctly classified (TP+TN)/(TP+TN+FP+FN). Finally, the Kappa Coefficient measures the agreement between the classification and the truth ((TN+FP) (TN+FN) + (FN+TP) (FP+TP)) / (TP+FP+TN+FN)2
}
\examples{

# Get hoopoe data and crop where appropriate
data("hoopoe")
start = as.POSIXct("2016-08-01","\%Y-\%m-\%d", tz="UTC")
end = as.POSIXct("2017-05-15","\%Y-\%m-\%d", tz="UTC")
PAM_data = create_crop(hoopoe, start, end)

# Generate a classification from activity
activity_classification = classify_flap(PAM_data$acceleration, period = 10)
reference = activity_classification$classification
reference[reference != activity_classification$migration] = "Not Migrating"
reference[(activity_classification$classification ==
          activity_classification$high_movement)] = "Flying"
reference[reference == activity_classification$migration] = "Migrating"

# Generate a very simple pressure classification based on pressure changes
classified  = ifelse( abs(diff(PAM_data$pressure$obs)) >2, "Migrating", "Not Migrating")

# Remove dates from activity which are not also present in pressure
to_keep = which(PAM_data$acceleration$date \%in\% PAM_data$pressure$date)
reference = reference[to_keep]
reference = reference[-1]

compare_confusion_matrix(reference, classified)

}
\references{
Congalton, R.G. and Green, K., 2008. Assessing the accuracy of remotely sensed data: principles and practices. CRC press.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/importPAM.R
\name{importPAM}
\alias{importPAM}
\title{Import PAM data}
\usage{
importPAM(
  pathname = pathname,
  measurements = c(".pressure", ".glf", ".acceleration", ".temperature", ".magnetic")
)
}
\arguments{
\item{pathname}{path where files are stored}

\item{measurements}{a series of measurements logged by the PAM logger which are to be imported. Currently supports these file extentions: ".pressure", ".glf", ".gle",".acceleration", ".temperature", "AirTemperature", ".BodyTemperature" and ".magnetic"}
}
\value{
a list of measurements for the one individual
}
\description{
Imports and formats many datasets into one big nested list containing all the data from the different sensors. A subset of sensors can be selected using \code{measurements}.
}
\examples{
#pathname = "your/filepath"
#measurements = c(".pressure", ".glf")
#PAM_data = importPAM(pathname, measurements)
#str(PAM_data)
#plot(PAM_data$light$date[3000:5000], PAM_data$light$obs[3000:5000],
#type="l", xlab="Date", ylab="Light Intensity")

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/create_merged_classification.R
\name{create_merged_classification}
\alias{create_merged_classification}
\title{Convert data from a Timetable to a timeseries}
\usage{
create_merged_classification(from, to, classification, add_to, missing = NA)
}
\arguments{
\item{from}{start of event that was classified (generally timetable output)}

\item{to}{end of event that was classified (generally timetable output)}

\item{classification}{classified data}

\item{add_to}{data which the classifications are to be added to (e.g. PAM_data$pressure)}

\item{missing}{Missing value replacement. By default NA.}
}
\value{
the classification in add_to dataset
}
\description{
convert a classification timetable into a classification timeseries
}
\examples{
data(bee_eater)
PAM_data = bee_eater

twl = GeoLight::twilightCalc(PAM_data$light$date, PAM_data$light$obs,
                             LightThreshold = 2, ask = FALSE)


to_classify = create_summary_statistics(PAM_data,
                     method="pressure",
                     twl = twl,
                     Pdiff_thld = 2,
                     light_thld = 2)

classification = classify_summary_statistics((to_classify$total_daily_duration *
                              log(to_classify$night_P_diff+0.001 )
                              * to_classify$total_daily_P_change),
                             states=3, "hmm")$cluster

pressure_classification = create_merged_classification(from = to_classify$start,
                                                       to =to_classify$end,
                                                       classification = classification,
                                                       add_to = PAM_data$pressure,
                                                       missing = NA)

pressure_classification[pressure_classification == NA] = 0

plot(PAM_data$pressure$date, PAM_data$pressure$obs,
     col= viridis::viridis(4)[pressure_classification+1],
     type="o", pch=16, cex=0.6)

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/confusionMAT.R
\name{confusionMAT}
\alias{confusionMAT}
\title{create a confusion matrix}
\usage{
confusionMAT(reference, classified)
}
\arguments{
\item{reference}{Reference dataset wich classification is to be compared to}

\item{classified}{Classification output to compare with a reference dataset}
}
\value{
a confusion matrix
}
\description{
The function populates a confusion matrix using predicted and reference points and estimate the errors in commission, omissioin for each class, in addition to Producer, User and overall accuracy, and the Kappa Coefficient. Errors in Commission provide a measure of false negatives i.e. the number of points that were predicted to be part of a class that they were not (probability something was incorrectly classified FN/(TP+FN)). Errors in Omission provide a measure of false positives that were predicted to be in a different class from their actual class (probability that something was missed FP/(FP +TP). User Accuracy or Recall represents the probability that a class was correctly classified TP/(TP + FN). Producer Accuracy or Precision provides a measure of how likely something was missed by the classification (probability that something was not missed TP/(TP + FP)). The Overall Accuracy represents the probability that all classes were correctly classified (TP+TN)/(TP+TN+FP+FN). Finally, the Kappa Coefficient measures the agreement between the classification and the truth ((TN+FP) (TN+FN) + (FN+TP) (FP+TP)) / (TP+FP+TN+FN)2
}
\examples{

## Get hoopoe data and crop where appropriate
#data("hoopoe")
#start = as.POSIXct("2016-08-01","\%Y-\%m-\%d", tz="UTC")
#end = as.POSIXct("2017-05-15","\%Y-\%m-\%d", tz="UTC")
#PAM_data = cutPAM(hoopoe, start, end)

## Generate a classification from activity
#activity_classification = classifyFLAP(PAM_data$acceleration, period = 10)
#reference = activity_classification$classification
#reference[reference != activity_classification$migration] = "Not Migrating"
#reference[(activity_classification$classification ==
#          activity_classification$high_movement)] = "Flying"
#reference[reference == activity_classification$migration] = "Migrating"

## Generate a very simple pressure classification based on pressure changes
#classified  = ifelse( abs(diff(PAM_data$pressure$obs)) >2, "Migrating", "Not Migrating")

## Remove dates from activity which are not also present in pressure
#to_keep = which(PAM_data$acceleration$date \%in\% PAM_data$pressure$date)
#reference = reference[to_keep]
#reference = reference[-1]

#confusionMAT(reference, classified)

}
\references{
Congalton, R.G. and Green, K., 2008. Assessing the accuracy of remotely sensed data: principles and practices. CRC press.
}
