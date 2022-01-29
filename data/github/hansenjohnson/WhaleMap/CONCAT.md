# WhaleMap
Collate and display whale survey results

[![DOI](https://joss.theoj.org/papers/10.21105/joss.03094/status.svg)](https://doi.org/10.21105/joss.03094)

## Overview

The goal of this software is to rapidly and effectively collect and share whale survey information within and between research, government, industry, and public sectors. Our hope is that it will improve survey efficiency, increase public awareness, and inform impactful, transparent management decisions. The system is live at [whalemap.org](https://whalemap.org/), and the interactive map is available [here](https://whalemap.org/WhaleMap). Check out the [WhaleMap publication in JOSS](https://doi.org/10.21105/joss.03094) for a general description of the system, or continue reading for more detailed information. Please report any problems or suggestions via this site's issue page. 

## Project structure

```
data/           All WhaleMap data (NOT tracked via git)
    raw/        Raw data cloned directly from data contributors (NEVER manually edited)
    interim/    Data from each platform coeherced to WhaleMap format
    processed/  Processed data used for display
R/              R scripts for data processing, display, and reporting
src/            Shell scripts to execute data cloning and processing
Makefile        Maps project dependency structure and facilitates efficient processing with `make`
LICENSE         License information
global.R        Shiny app
server.R        Shiny app
ui.R            Shiny app
```

## Data processing

### Workflow

1. Sync

The majority of the remoate data are synced using the script `src/get_remote_data.sh`. This uses Rclone to sync data from remote repositories (Google Drive, Dropbox, etc.), and also calls `src/get_live_dcs.sh` to download acoustic detection data. All the data are stored in `data/raw/`

2. Process 

Data from each contributor are processed to a common WhaleMap format (see below) using a custom R script (`R/proc_*.R`). Once formatted, observations and effort data from each platform are saved in `data/interim/`.

3. Combine 

All track data files in `data/interim/` are combined by (`R/proc_tracks.R`) and saved as `data/processed/tracks.rds`. All observation data files in `data/interim/` are combined by (`R/proc_observations.R`) and saved as `data/processed/observations.rds`.

4. Repeat

A makefile maps the dependency structure and orchestrates the efficient processing of the entire dataset. The `make` command is executed at the final step of `src/get_remote_data.sh` to update the dataset after synchronization. A cron job runs `src/get_remote_data.sh` every 15 minutes to keep the system up to date.

### WhaleMap data formats

#### Observations

`time` - UTC time (YYYY-MM-DD HH:MM:SS)  
`lat` - latitude in decimal degrees  
`lon` - longitude in decimal degrees  
`date` - UTC date (YYYY-MM-DD)  
`yday` - day of year  
`year` - year (YYYY)  
`platform` - type of survey platform (`vessel`, `plane`, `slocum`, `buoy`, `rpas`, `opportunistic`)  
`name` - name of platform (e.g., `noaa_twin_otter`, `dfo_coriolis`)  
`id` - unique survey identifier comprised of survey start date, platform, and name (e.g., `2020-02-21_plane_noaa_twin_otter`)  
`species` - species name (`right`, `fin`, `sei`, `humpback`, `blue`)  
`score` - detection type and score (`Definite acoustic`, `Possible acoustic`, `Definite visual`, `Possible visual`)  
`number` - number of whales (`NA` for acoustic detections)  
`calves` - number of calves (`NA` for acoustic detections)  

#### Tracks

`time` - UTC time (YYYY-MM-DD HH:MM:SS)  
`lat` - latitude in decimal degrees  
`lon` - longitude in decimal degrees  
`date` - UTC date (YYYY-MM-DD)  
`yday` - day of year  
`year` - year (YYYY)  
`platform` - type of survey platform (`vessel`, `plane`, `slocum`, `buoy`, `rpas`, `opportunistic`)  
`name` - name of platform (e.g., `noaa_twin_otter`, `dfo_coriolis`)  
`id` - unique survey identifier comprised of survey start date, platform, and name (e.g., `2020-02-21_plane_noaa_twin_otter`)  
`speed` - platform speed (m/s)  
`altitude` - platform altitude (m)  

#### Status table

`script` - name of the platform processing script (e.g., `proc_2021_noaa_twin_otter.R`)  
`name` - name of the platform to be displayed in the status table (e.g., `NOAA NEFSC survey sightings/tracks`)  
`url` - link url to be displayed in status table  
`email_list` - path to csv file with list of emails to be notified if there is an error processing platform data  
 
### Reporting

Error catching is performed in the `Makefile` using the scripts `src/report_error.sh` and `src/remove_error.sh`. These are run each time a platform-specific processing script is run. The results of the processing are recorded in the status table (`status.txt`). If processing is successful, a timestamp is added in the status table. If processing is unsuccessful, an error message is printed in the status table and an auto-generated email is sent to a designated email list associated with that platform (using `src/send_email_alert.sh`). The status table is displayed in the Shiny app.

### Adding a new platform

1. Update `src/get_remote_data.sh` to sync raw data  
2. Write R script to convert to WhaleMap format and save observations and tracks in `data/interim.`  
3. Update `Makefile`  
4. Update `status.txt`  
5. Test (run `make` and view results in Shiny app)  
6. Push changes to WhaleMap server via GitHub  

## Setup

### R

Here is a list of the packages that WhaleMap relies on:
- Data wrangling: `tidyverse`, `lubridate`, `tools`, `oce`   
- Mapping: `rgdal`, `maptools`, `leaflet`, `leaflet.extras`, `sf`   
- Shiny: `shiny`, `shinydashboard`, `shinybusy`   
- Misc: `htmltools`, `htmlwidgets`, `plotly`, `RColorBrewer`

Here's the output from `sessionInfo()`:
```
> sessionInfo()
R version 3.4.4 (2018-03-15)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: Ubuntu 16.04.6 LTS

Matrix products: default
BLAS: /usr/lib/libblas/libblas.so.3.6.0
LAPACK: /usr/lib/lapack/liblapack.so.3.6.0

locale:
 [1] LC_CTYPE=en_CA.UTF-8       LC_NUMERIC=C              
 [3] LC_TIME=en_CA.UTF-8        LC_COLLATE=en_CA.UTF-8    
 [5] LC_MONETARY=en_CA.UTF-8    LC_MESSAGES=en_CA.UTF-8   
 [7] LC_PAPER=en_CA.UTF-8       LC_NAME=C                 
 [9] LC_ADDRESS=C               LC_TELEPHONE=C            
[11] LC_MEASUREMENT=en_CA.UTF-8 LC_IDENTIFICATION=C       

attached base packages:
[1] tools     stats     graphics  grDevices utils     datasets  methods  
[8] base     

other attached packages:
 [1] rgeos_0.3-26         RColorBrewer_1.1-2   forcats_0.4.0       
 [4] stringr_1.3.1        dplyr_0.8.0.1        purrr_0.3.0         
 [7] readr_1.1.1          tidyr_0.8.2          tibble_2.0.1        
[10] tidyverse_1.2.1      sf_0.7-3             shinybusy_0.1.2     
[13] leaflet.extras_1.0.0 plotly_4.8.0         ggplot2_3.1.0       
[16] shinydashboard_0.7.0 oce_0.9-23           gsw_1.0-5           
[19] testthat_2.0.0       lubridate_1.7.4      maptools_0.9-2      
[22] htmlwidgets_1.3      htmltools_0.3.6      rgdal_1.2-18        
[25] sp_1.3-1             leaflet_2.0.1        shiny_1.2.0  
```

### Linux / tools

Here's a list of other tools that WhaleMap relies on and the current version:
- make (GNU Make 4.1)
- rclone (rclone v1.38)
- git (git version 2.7.4)
- pandoc (pandoc 1.16.0.2)
- shiny-verser (Shiny Server v1.5.7.907)# Syntax examples
***

## chmod
Make a shell script executable

```
# add the following to the first line of the file
#!/bin/bash

# run this line to make executable
chmod u+x scriptname.sh
```

***

## rsync
File sync / transfer

```
# transfer to WhaleMap
rsync -rtv file.txt hansen@whalemapvm:/srv/shiny-server/WhaleMap/
```
***

## crontab
Automate tasks

```
# open crontab
crontab -e

# enable edits to crontab
a # this means 'append' in vim

# run every hour
0 * * * * sh /home/hansen/shiny-server/WhaleMap/get_live_dcs.sh

# run every 5 min (usually just for testing)
*/5 * * * * sh /home/hansen/shiny-server/WhaleMap/get_live_dcs.sh

# server side example
0 * * * * /srv/shiny-server/WhaleMap/get_remote_data.sh

# quit and save crontab
:x
```
***
## Shiny server

```
# reload the shiny server
sudo systemctl reload shiny-server
# restart the shiny server
sudo systemctl restart shiny-server
```

***

## Cron mail

```
# Delete
cat /dev/null >/var/mail/hansenjohnson # local
cat /dev/null >/var/mail/hansen # server
```
***

## Git

### Delete local commits
```
git reset --hard
```

### Apply patch commit
```
git am filename
```

#### Merge branches (dev into master)
https://git-scm.com/book/en/v2/Git-Branching-Basic-Branching-and-Merging
```
git checkout master
git merge dev
```
---
title: 'WhaleMap: a tool to collate and display whale survey results in near real-time'
tags:
  - R
  - oceanography
  - marine mammals
  - conservation
authors:
  - name: Hansen Johnson
    orcid: 0000-0002-3086-6759
    affiliation: 1
  - name: Daniel Morrison
    affiliation: 1
  - name: Christopher Taggart
    affiliation: 1
affiliations:
 - name: Oceanography Department, Dalhousie University, 1355 Oxford Street, Halifax, Nova Scotia, Canada B3H 4R2
   index: 1
date: 18 February 2021
bibliography: paper.bib
---

# Statement of Need

Baleen whales of the Northwest Atlantic live in a highly urbanized ocean. Their recovery from commercial whaling is impeded by anthropogenic risks from ocean industry, pollution, and climate change. Effective research, conservation and risk-reduction action requires near real-time knowledge of whale distribution measured using various methods including visual surveys from vessels or planes or acoustic surveys from autonomous platforms. The rapid collation and dissemination of whale detections and survey effort is critical but challenging given the number and variety of survey organizations and methodologies at work along the east coast of the US and Canada. There are long term databases for whale survey data, such as that maintained by the North Atlantic Right Whale Consortium ([narwc.org](https://narwc.org)), and crowd-source reporting tools (e.g., [Whale Alert](http://www.whalealert.org/)) but `WhaleMap` is the only dedicated system specifically designed to collate and display all available near real-time whale detections and survey effort. Use cases vary widely. For example, `WhaleMap` is currently used by: government managers to design and implement risk-mitigation strategies, members of military or industry to plan safe operations, researchers to coordinate survey efforts and explore patterns in whale distribution, and members of the general public to learn about and follow along with whale conservation activities.

`WhaleMap` was designed with several specific goals:  
-	Incorporate whale detection and survey effort from all survey methods in near real-time  
-	Allow survey teams to easily contribute and retain complete control over their data  
-	Provide the latest data in an accurate, user-friendly, and publicly accessible format  
-	Operate transparently using open-source tools and with limited supervision  

Critically, `WhaleMap` does **not**:  
-	Perform any quality-control, or take responsibility for the veracity of information contributed  
-	Provide a long-term database for survey results  
-	Allow access to raw or processed data without approval from the data originator  

# System

The `WhaleMap` system workflow can be separated into data processing and visualization components (Figure 1). The following provides a brief overview of each. Additional details as well as specific references to all software used is available in the source code documentation.

## Data processing

Survey teams provide `WhaleMap` access to a remote repository of their choice (e.g., Google Drive, Dropbox) where they upload their survey data. The `WhaleMap` curator writes a custom script to extract the detection and effort data from each survey team and convert it to a common `WhaleMap` format. This method eases the burden on the survey teams by allowing any team to submit data in nearly any format, provided the format is consistent and well-documented. This is essential for rapid data collection, as survey teams in the field typically lack the time and resources to reformat their data.

A scheduled job regularly clones the data from the remote repositories onto the `WhaleMap` server and uses a makefile to dynamically and efficiently process the data from each platform and coerce it into a common format. Formatting errors in a remote data repository are automatically flagged and the contributor is notified. This ensures that any changes to the raw survey data quickly propagate through the entire system, which allows survey teams to retain complete control of their data and perform quality control as needed. It also guarantees that the `WhaleMap` always contains the latest available information.

## Visualization

Once the survey data are processed, they are visualized using two different methods. The first is the construction of self-contained HTML summary maps containing sufficient information to satisfy most casual viewers (typically the last 14-days of survey results). These can be conveniently embedded in various webpages (e.g., [whalemap.org](https://whalemap.org)) and browsed without requiring server-side processing. These maps are dynamically regenerated as the final step in the data processing workflow, so they always contain the latest available information. The second visualization method is an interactive online application ([whalemap.org/WhaleMap](https://whalemap.org/WhaleMap/)). This provides users with numerous tools with which to filter the latest processed data. The selected data are displayed in several formats including an interactive map, interactive timeseries plot, and table of summary statistics. 

# Conclusions

Since its launch in 2018, `WhaleMap` has been constantly refined and optimized to better serve the overall goal of providing a common source for all near real-time whale survey data in the Northwest Atlantic. It has demonstrably improved conservation outcomes for endangered whales in this region by optimizing research activities, facilitating dynamic risk-mitigation measures, and engaging with the ocean industry and the public. `WhaleMap` has also already been cited in several scientific publications [@gervaise:2021; @koubrak:2021; @baumgartner:2020; @johnson:2020; @kowarski:2020]. We are not aware of any equivalent software in existence. It is our hope that `WhaleMap` continues to serve the conservation community in perpetuity, and that this system can be readily adapted to benefit other regions facing similar conservation challenges.

![Conceptual overview of WhaleMap data processing and display](figure_1.png)

# Acknowledgements

We gratefully acknowledge the support of Pamela Emery, Stephanie Ratelle, Angelia Vanderlaan, Hilary Moors-Murphy and many other colleagues at Fisheries and Oceans Canada, as well as Tim Cole, Elizabeth Josephson, Leah Crowe, Christin Khan, Danielle Cholewiak and others at the NOAA Northeast Fisheries Science Center. This work would not have been possible without them. We also thank Mark Baumgartner, Kim Davies, and members of the Taggart lab for advice and helpful conversations. We are indebted to the fantastic community of open source developers that constructed many of the tools on which `WhaleMap` relies. These are referenced in detail in the system documentation. Finally, we give our most sincere thanks to the many survey teams and numerous agencies that continue to place their trust in this system and praise their tireless efforts to protect these vulnerable species. Supporting funds generously provided by Fisheries and Oceans Canada and the Natural Sciences and Engineering Research Council of Canada. 

# References
