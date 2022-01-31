# phonfieldwork

[![Project Status: Active - The project has reached a stable, usable state and is being actively developed.](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active)
[![](https://badges.ropensci.org/385_status.svg)](https://github.com/ropensci/software-review/issues/385)
[![CRAN version](http://www.r-pkg.org/badges/version/phonfieldwork)](https://cran.r-project.org/package=phonfieldwork)
[![](http://cranlogs.r-pkg.org/badges/grand-total/phonfieldwork)](https://CRAN.R-project.org/package=phonfieldwork)
[![R build status](https://github.com/ropensci/phonfieldwork/workflows/R-CMD-check/badge.svg)](https://github.com/ropensci/phonfieldwork/actions)
[![Coverage Status](https://img.shields.io/codecov/c/github/ropensci/phonfieldwork/master.svg)](https://codecov.io/github/ropensci/phonfieldwork?branch=master)
[![DOI](https://zenodo.org/badge/194053227.svg)](https://zenodo.org/badge/latestdoi/194053227)

`phonfieldwork` is a package for phonetic fieldwork research and experiments. This package make it easier:

- creating a html/pptx presentation from stimuli-translation list, 
- renaming soundfiles according to the list of stimuli, 
- concatenating multiple soundfiles and create a Praat TextGrid whose interval labels are the original names of the sound
- extracting sounds according to annotation
- extracting annotation from multiple linguistic formats (Praat `.TextGrid`, ELAN `.eaf`, EXMARaLDA `.exb`, Audacity `.txt` and subtitles `.srt`)
- visualising an oscilogram, a spectrogram and an annotation
- creating an html viewer [like this](https://ropensci.github.io/phonfieldwork/s1/stimuli_viewer.html), ethical problems of this kind of viewer in linguistic research are covered in the vignette `vignette("ethical_research_with_phonfieldwork")`.

For more ditails see [tutorial](https://docs.ropensci.org/phonfieldwork/).

The main goal of the `phonfieldwork` package is to make a full research workflow from data collection to data extraction and data representation easier for people that are not familiar with programming. Hovewer most of the `phonfieldwork` funnctionality can be found in other software and packages:

* stimuli presentation creation could be done with any programming language and probably without them
* automatic file renaming and automatic merge could be done with any programming language
* Praat `.TextGrid` manipulation could be done with Praat, R packages [`rPraat`](https://cran.r-project.org/package=rPraat) and [`textgRid`](https://cran.r-project.org/package=textgRid), the Python package ['pympi'](https://dopefishh.github.io/pympi/index.html))
* ELAN `.eaf` manipulationcould be done with ELAN, the R package [`FRelan`](https://github.com/langdoc/FRelan) and the Python package [`pympi`](https://dopefishh.github.io/pympi/index.html)
* Praat `.TextGrid`, ELAN `.eaf`, and 'EXMARaLDA .exb import and export could be done with the R package [`act`](https://cran.r-project.org/package=act)
* cut sounds according to annotation could be done with Praat and the R package`tuneR`
* spectrogram visualisation could be done with multiple R packages [`signal`](https://cran.r-project.org/package=signal), [`tuneR`](https://cran.r-project.org/package=tuneR), [`seewave`](https://cran.r-project.org/package=seewave), [`phonTools`](https://cran.r-project.org/package=phonTools), [`monitor`](https://cran.r-project.org/package=monitor), [`warbleR`](https://cran.r-project.org/package=warbleR), [`soundgen`](https://cran.r-project.org/package=soundgen) and many others

## Installation

Install from CRAN:

```
install.packages("phonfieldwork")
```

Get the development version from GitHub:

```
install.packages("remotes")
remotes::install_github("ropensci/phonfieldwork")
```
Load a library:
```
library(phonfieldwork)
```

In order to work with some `rmarkdown` functions you will need to install `pandoc`, see `vignette("pandoc")` for the details.

## To do:

* export to ELAN and EXMARALDA files
* use ELAN and EXMARALDA files in the whole pipline discribed in docs
* use the same pipline with video (for Sign Languages)
* make [TECkit](https://scripts.sil.org/cms/scripts/render_download.php?format=file&media_id=BeyondUTR22_pdf&filename=BeyondUTR22_pdf.pdf) to df and back

[![ropensci\_footer](https://ropensci.org/public_images/ropensci_footer.png)](https://ropensci.org)
# phonfieldwork 0.0.12

- add a `separate_duration` argument to the `concatenate_soundfiles()` function that makes it possible to use some silent separator during the file concatenation.
- make `rename_soundfiles()` function to work with mp3 files.
- fix the bug with `"` sign in textgrids.
- make pictures optional in the `create_viewer()` function

# phonfieldwork 0.0.11

- correct empty tiers behavior #34 (thanks to Shungo Suzuki)
- add possibility to have different values in the `n_of_annotations` argument of `create_subannotation()` (thanks to Jenya Korovina for the idea)
- rename `tier` argument of the `create_empty_textgrid()` to `tier_name`.
- create the `remove_textgrid_tier()` function.

# phonfieldwork 0.0.10

- add tryCatch to the `read_from_folder()` function

# phonfieldwork 0.0.8

- add possibility to read short format of `.TextGrid`s and fix `textgrid_to_df()` and `tier_to_df()` functions
- add `encoding`, `formant_df`, `intensity`, `picth` and `pitch_range` arguments to the `draw_sound()` function
- add the `formant_to_df()` function
- add the `picth_to_df()` function
- add the `intensity_to_df()` function
- add an argument `external` to the `create_presentation()` function in order to mark external images or gifs
- remove all `encoding` arguments and replace it with encoding autodetection from `uchardet` (thanks to Artem Klevtsov for help)
- add `autonumber`, `loging` and `missing` arguments to `rename_soundfiles()` function (thanks to Niko Partanen)
- a lot of minor style changes (thanks to Jonathan Keane)
- add the `create_empty_textgrid()` function (thanks to Niko Partanen)
- add the `data_manipulation_with_tidyverse` vignette (thanks to Niko Partanen)
- add the `concatenate_textgrids()` function
- add the `read_from_folder()` function and remove `..._from_folder` arguments
- pass rOpenSci review! Move tutorial to <https://ropensci.github.io/phonfieldwork/>

# phonfieldwork 0.0.7

- add a vigniettes about ethical research and introduction to work with phonfieldwork
- add an argument `textgrids_from_folder` to the `textgrid_to_df()` function
- add an argument `exbs_from_folder` to the `exb_to_df()` function
- add an argument `eafs_from_folder` to the `eaf_to_df()` function
- add Raven style annotations8
- replace freqmax with frequency_range argument
- replace example_textgrid with systemfile() call
- add window annotation to spectrograms
- add bridge to lingtypology package: `map` argument in the `create_viewer()` function
- make it possible to visualise all types of annotations with the `draw_sound()` function
- add the `source` column to all `..._to_df()` functions
- add the `audacity_to_df` function
- add the `srt_to_df` function
- chage textgrid related functions' output from `start`, `end`, `annotation` to `time_start`, `time_end`, `content`
- correct point tier visualization

# phonfieldwork 0.0.6

- add `encoding` arguments to functions for working with TextGrids
- add `textgrid` argument to the `draw_sound()` function
- change subgraphs alignment in the `draw_sound()` function including textgrid annotation
- add `from` and `to` arguments to the `draw_sound()` function
- add .mp3 format reading options to all functions that work with sounds
- add `text_size` argument to the `draw_sound()` function
- add `zoom`` argument to the `draw_sound()` function
- add the `get_sound_duration()` function
- fix ploting of multiple sounds with multiple .TextGrids
- add an argument `title_as_filename` to the `draw_sound()` function
- change .TextGrid associated arguments of the `create_viewer()` function to `table` argument; as a result users now need to provide a table for the annotation viewer and not a .TextGrid

# phonfieldwork 0.0.5

- add `textgrid_to_df()` function for reading Praat files
- add `create_glossed_document()` function for converting .flextext files into a glossed document
- add `flextext_to_df()` function for reading FLEx files
- add `eaf_to_df()` function for reading ELAN files
- add `exb_to_df()` function for reading EXMARaLDA files
- rename `textgrid` argument into `annotation` argument in `concatenate_soundfiles()` function adding new possible values

# phonfieldwork 0.0.4

- add `create_subannotation()` function

# phonfieldwork 0.0.3

- vertically and horisontally center text in presentations created by `create_presentation()`; thx @Pandaklez #1
- add the `font_size` argument to the `create_presentation()` function
- add `rename_videofiles()` function
- rebuild html viewer for sounds with JavaScript with the help of new functions `create_image_look_up()` and `create_sound_play()`.

# phonfieldwork 0.0.2

- make the `create_presentation()` function render silently
- add a new function `draw_sound()` for creating spectrogram and oscilogram
- add a new function `create_viewer()` for creating an html viewer with sound and spectrograms
- correct work of `autonumbering` function in `extract_intervals()` function
- add new functions  `get_textgrid_names()` and `set_textgrid_names()`
- finish tutorial <https://ropensci.github.io/phonfieldwork/>

# phonfieldwork 0.0.1

- initial release
# Contributing to `phonfieldwork`

## Issues

When filing an issue, the most important thing is to include a minimal reproducible example so that I can quickly verify the problem. So please include:

* required packages
* package versions
* data
* code

There are some additional information that I put in the issue template.
---
name: Bug report
about: Create a report to help us improve
title: ''
labels: bug
assignees: ''
---

Please, put the results of these functions below, if your issue is related to the technical bugs. Thank you!

<details> <summary> info about OS and package versions </summary>
```
sessionInfo()$R.version$platform
sessionInfo()$R.version$version.string
packageVersion("rmarkdown")
packageVersion("phonfieldork")
```
</details>
---
name: Feature request
about: 'Suggest an idea for this project '
title: ''
labels: enhancement
assignees: ''

---


---
title: "Sound Viewer"
output: 
  html_document:
      self_contained: yes
params:
  data: your_file.csv
  captions: "caption_for_pictures"
  about: about.Rmd
  map: map.Rmd
editor_options: 
  chunk_output_type: console
---

<style>
#my_block {
  display: none;
  position: fixed;
  z-index: 1;
  padding-top: 100px; 
  left: 0;
  top: 0;
  width: 100%; 
  height: 100%;
  overflow: auto; 
  background-color: rgba(0,0,0, 0.1);
}
#my_img {
  margin: auto;
  display: block;
  width: 80%;
  max-width: 700px;
}
#caption {
  margin: auto;
  display: block;
  width: 80%;
  max-width: 700px;
  text-align: center;
  padding: 10px 0;
  height: 150px;
}
</style>

```{r data_creation, include=FALSE}
df <- read.csv(params$data, stringsAsFactors = FALSE)

df$viewer <- 
  paste(create_image_look_up(img_src = df$pictures,
                                            img_caption = params$captions),
        create_sound_play(snd_src = df$audio))

df <- df[, -which(colnames(df) %in% c("audio", "pictures"))]
```


# {.tabset .tabset-fade .tabset-pills}

```{r add_map, child = params$map}
```

## data

```{r make_viewer, echo=FALSE, message=FALSE, warning=FALSE}
DT::datatable(data = df, 
              filter = 'top', 
              rownames = FALSE, 
              options = list(pageLength = 50, dom = 'tip'), 
              escape = FALSE)
```

<div class = "my_block" id="my_block" onclick = "pic_disappear()">
  <span class="close">&times;</span>
  <img class = "my_img" id="my_img">
  <div class = "caption" id="caption">
  </div>
</div>


## about

```{r about, child = params$about}
```


<script>
  function pic_appear(path, caption) {
    var path, caption;
    var block = document.getElementById("my_block");
    var block_img = document.getElementById("my_img");
    var caption_text = document.getElementById("caption");
    block.style.display = "block";
    block_img.src = path;
    caption_text.innerHTML = caption;
  }
  function pic_disappear() {
    var block = document.getElementById("my_block");
    block.style.display = "none";
  }
  function sound_play(x) {
    var audio = new Audio();
    audio.src = x;
    audio.play();
  }
  function resize(elem, percent) { 
  elem.style.fontSize = percent; 
  }
</script>
---
title: "Glossed document by phonfieldwork [Moroz 2020]"
params:
  data: your_file.csv
  rows: rows
  example_pkg: example_pkg
editor_options: 
  chunk_output_type: console
---

```{r make_viewer, echo=FALSE, message=FALSE, warning=FALSE, results='asis'}
df <- read.csv(params$data, stringsAsFactors = FALSE)

all_rows <- c("cf", "hn", "gls", "msa")

library(dplyr, warn.conflicts = FALSE)
library(tidyr, warn.conflicts = FALSE)

if(is.null(example_pkg)){
  df %>%
    group_by(text_title, p_id, s_id, free_trans, w_id) %>%
    summarise(txt = paste0(txt, collapse = ""),
              cf = paste0(cf, collapse = "-"),
              hn = paste0(hn, collapse = "-"),
              msa = paste0(msa, collapse = "-"),
              gls = paste0(gls, collapse = "-")) %>% 
    mutate(cf = gsub("-=", "=", cf),
           hn = gsub("-=", "=", hn),
           msa = gsub("-=", "=", msa),
           gls = gsub("-=", "=", gls)) %>% 
    group_by(text_title, p_id, s_id, free_trans) %>%
    summarise(txt = paste0(txt, collapse = "&#x09;"),
              cf = paste0(cf, collapse = "&#x09;"),
              hn = paste0(hn, collapse = "&#x09;"),
              msa = paste0(msa, collapse = "&#x09;"),
              gls = paste0(gls, collapse = "&#x09;")) %>% 
    mutate(txt = paste0(s_id, "\\.&#x09;", txt, "\n\n&#x09;"),
           cf = paste0(cf, "\n\n&#x09;"),
           hn = paste0(hn, "\n\n&#x09;"),
           msa = paste0(msa, "\n\n&#x09;"),
           gls = paste0(gls, "\n\n&#x09;"),
           free_trans = paste0(free_trans, "\n\n")) %>% 
    select(text_title, p_id, s_id, txt, cf, hn, msa, gls, free_trans) %>% 
    ungroup() %>% 
    add_row(tibble(text_title = unique(.$text_title),
                   p_id = 0,
                   s_id = 0,
                   txt = paste0("\n\n## ", unique(.$text_title), "\n\n"),
                   cf = NA,
                   hn = NA,
                   msa = NA,
                   gls = NA,
                   free_trans = NA), .before = 1)  %>% 
    arrange(text_title, p_id, s_id) %>% 
    select(-c(all_rows[!(all_rows %in% params$rows)])) %>% 
    pivot_longer(names_to = "type", values_to = "value", cols = c("txt", params$rows, "free_trans")) %>% 
    na.omit() %>% 
    ungroup() %>% 
    select(value) %>% 
    unlist() %>% 
    unname() %>% 
    cat(sep = "")
} else if(example_pkg == "philex"){
  df %>%
    group_by(text_title, p_id, s_id, free_trans, w_id) %>%
    summarise(txt = paste0(txt, collapse = ""),
              cf = paste0(cf, collapse = "-"),
              hn = paste0(hn, collapse = "-"),
              msa = paste0(msa, collapse = "-"),
              gls = paste0(gls, collapse = "-")) %>% 
    mutate(cf = gsub("-=", "=", cf),
           hn = gsub("-=", "=", hn),
           msa = gsub("-=", "=", msa),
           gls = gsub("-=", "=", gls)) %>%
    group_by(text_title, p_id, s_id, free_trans) %>%
    summarise(txt = paste0(txt, collapse = " "),
              cf = paste0(cf, collapse = " "),
              hn = paste0(hn, collapse = " "),
              msa = paste0(msa, collapse = " "),
              gls = paste0(gls, collapse = " ")) %>%
    mutate(txt = paste0("\\\\lb\\{ex:", s_id, "\\}\\{\\\\gll ", txt, "\\\\\\\\\n\n&#x09;"),
           cf = paste0(cf, "\\\\\\\\\n\n&#x09;"),
           hn = paste0(hn, "\\\\\\\\\n\n&#x09;"),
           msa = paste0(msa, "\\\\\\\\\n\n&#x09;"),
           gls = paste0(gls, "\\\\\\\\\n\n&#x09;"),
           free_trans = paste0("\\\\trans `", free_trans, "'\\}\n\n")) %>% 
    select(text_title, p_id, s_id, txt, cf, hn, msa, gls, free_trans) %>% 
    ungroup() %>% 
    add_row(tibble(text_title = unique(.$text_title),
                   p_id = 0,
                   s_id = 0,
                   txt = paste0("\n\n## ", unique(.$text_title), "\n\n"),
                   cf = NA,
                   hn = NA,
                   msa = NA,
                   gls = NA,
                   free_trans = NA), .before = 1)  %>% 
    arrange(text_title, p_id, s_id) %>% 
    select(-c(all_rows[!(all_rows %in% params$rows)])) %>%
    pivot_longer(names_to = "type", values_to = "value", cols = c("txt", params$rows, "free_trans")) %>%
    na.omit() %>% 
    ungroup() %>% 
    select(value) %>% 
    unlist() %>% 
    unname() %>% 
    cat(sep = "")
} else if(example_pkg == "gb4e"){
  # \begin{exe}
  # \ex\label{ex1}
  # \gll Den Fritz_1 habe ich \_\_{}_1 zum Essen eingeladen.\\
  # the fred have I {} {to the} eating invited.\\
  # \glt I invited Fred for dinner.
  # \end{exe}
  
  df %>%
    group_by(text_title, p_id, s_id, free_trans, w_id) %>%
    summarise(txt = paste0(txt, collapse = ""),
              cf = paste0(cf, collapse = "-"),
              hn = paste0(hn, collapse = "-"),
              msa = paste0(msa, collapse = "-"),
              gls = paste0(gls, collapse = "-")) %>% 
    mutate(cf = gsub("-=", "=", cf),
           hn = gsub("-=", "=", hn),
           msa = gsub("-=", "=", msa),
           gls = gsub("-=", "=", gls)) %>% 
    group_by(text_title, p_id, s_id, free_trans) %>%
    summarise(txt = paste0(txt, collapse = " "),
              cf = paste0(cf, collapse = " "),
              hn = paste0(hn, collapse = " "),
              msa = paste0(msa, collapse = " "),
              gls = paste0(gls, collapse = " ")) %>% 
    mutate(txt = paste0("\\\\begin\\{exe\\}\n\n",
                        "\\\\ex\\\\label\\{ex:", s_id, "\\}\n\n",
                        "\\\\gll ", txt, "\\\\\\\\\n\n&#x09;"),
           cf = paste0(cf, "\\\\\\\\\n\n&#x09;"),
           hn = paste0(hn, "\\\\\\\\\n\n&#x09;"),
           msa = paste0(msa, "\\\\\\\\\n\n&#x09;"),
           gls = paste0(gls, "\\\\\\\\\n\n&#x09;"),
           free_trans = paste0("\\\\glt `", free_trans, "'\n\n",
                               "\\\\end\\{exe\\}\n\n")) %>% 
    select(text_title, p_id, s_id, txt, cf, hn, msa, gls, free_trans) %>% 
    ungroup() %>% 
    add_row(tibble(text_title = unique(.$text_title),
                   p_id = 0,
                   s_id = 0,
                   txt = paste0("\n\n## ", unique(.$text_title), "\n\n"),
                   cf = NA,
                   hn = NA,
                   msa = NA,
                   gls = NA,
                   free_trans = NA), .before = 1)  %>% 
    arrange(text_title, p_id, s_id) %>% 
    select(-c(all_rows[!(all_rows %in% params$rows)])) %>% 
    pivot_longer(names_to = "type", values_to = "value", cols = c("txt", params$rows, "free_trans")) %>% 
    na.omit() %>% 
    ungroup() %>% 
    select(value) %>% 
    unlist() %>% 
    unname() %>% 
    cat(sep = "")
} else if(example_pkg == "langsci"){
  # \ea\label{ex:examplelabel}
  # \langinfo{French}{Indo-European}{personal knowledge}\\
  # \gll Jean aim-e
  # Marie \\
  # John love-\textsc{3s.pres.ind} Mary \\
  # \glt ‘John loves Mary.’
  # \z
  
  df %>%
    group_by(text_title, p_id, s_id, free_trans, w_id) %>%
    summarise(txt = paste0(txt, collapse = ""),
              cf = paste0(cf, collapse = "-"),
              hn = paste0(hn, collapse = "-"),
              msa = paste0(msa, collapse = "-"),
              gls = paste0(gls, collapse = "-")) %>% 
    mutate(cf = gsub("-=", "=", cf),
           hn = gsub("-=", "=", hn),
           msa = gsub("-=", "=", msa),
           gls = gsub("-=", "=", gls)) %>% 
    group_by(text_title, p_id, s_id, free_trans) %>%
    summarise(txt = paste0(txt, collapse = " "),
              cf = paste0(cf, collapse = " "),
              hn = paste0(hn, collapse = " "),
              msa = paste0(msa, collapse = " "),
              gls = paste0(gls, collapse = " ")) %>% 
    mutate(txt = paste0("\\\\ea\\\\label\\{ex:", s_id, "\\}\n\n",
                        "\\\\gll ", txt, "\\\\\\\\\n\n&#x09;"),
           cf = paste0(cf, "\\\\\\\\\n\n&#x09;"),
           hn = paste0(hn, "\\\\\\\\\n\n&#x09;"),
           msa = paste0(msa, "\\\\\\\\\n\n&#x09;"),
           gls = paste0(gls, "\\\\\\\\\n\n&#x09;"),
           free_trans = paste0("\\\\glt `", free_trans, "'\n\n",
                               "\\\\z\n\n")) %>% 
    select(text_title, p_id, s_id, txt, cf, hn, msa, gls, free_trans) %>% 
    ungroup() %>% 
    add_row(tibble(text_title = unique(.$text_title),
                   p_id = 0,
                   s_id = 0,
                   txt = paste0("\n\n## ", unique(.$text_title), "\n\n"),
                   cf = NA,
                   hn = NA,
                   msa = NA,
                   gls = NA,
                   free_trans = NA), .before = 1)  %>% 
    arrange(text_title, p_id, s_id) %>% 
    select(-c(all_rows[!(all_rows %in% params$rows)])) %>% 
    pivot_longer(names_to = "type", values_to = "value", cols = c("txt", params$rows, "free_trans")) %>% 
    na.omit() %>% 
    ungroup() %>% 
    select(value) %>% 
    unlist() %>% 
    unname() %>% 
    cat(sep = "")  
  
} else if(example_pkg == "expex"){
  # \ex\label{ex:exmaple1}\begingl
  # \gla k- wapm -a -s’i -m -wapunin -uk //
  # \glb CL V AGR NEG AGR TNS AGR //
  # \glb 2 see {\sc 3acc} {} {\sc 2pl} preterit {\sc 3pl} //
  # \glft ‘you (pl) didn’t see them’//
  # \endgl  \xe
  
  df %>%
    group_by(text_title, p_id, s_id, free_trans, w_id) %>%
    summarise(txt = paste0(txt, collapse = ""),
              cf = paste0(cf, collapse = "-"),
              hn = paste0(hn, collapse = "-"),
              msa = paste0(msa, collapse = "-"),
              gls = paste0(gls, collapse = "-")) %>% 
    mutate(cf = gsub("-=", "=", cf),
           hn = gsub("-=", "=", hn),
           msa = gsub("-=", "=", msa),
           gls = gsub("-=", "=", gls)) %>% 
    group_by(text_title, p_id, s_id, free_trans) %>%
    summarise(txt = paste0(txt, collapse = " "),
              cf = paste0(cf, collapse = " "),
              hn = paste0(hn, collapse = " "),
              msa = paste0(msa, collapse = " "),
              gls = paste0(gls, collapse = " ")) %>% 
    mutate(txt = paste0("\\\\ex\\\\label\\{ex:", s_id, "\\}//\begingl\n\n",
                        "\\\\gla ", txt, "//\n\n&#x09;"),
           cf = paste0("\\\\glb ", cf, "//\n\n&#x09;"),
           hn = paste0("\\\\glb ", hn, "//\n\n&#x09;"),
           msa = paste0("\\\\glb ", msa, "//\n\n&#x09;"),
           gls = paste0("\\\\glb ",gls, "//\n\n&#x09;"),
           free_trans = paste0("\\\\glft `", free_trans, "'//\n\n",
                               "\\\\endgl\\\\xe")) %>% 
    select(text_title, p_id, s_id, txt, cf, hn, msa, gls, free_trans) %>% 
    ungroup() %>% 
    add_row(tibble(text_title = unique(.$text_title),
                   p_id = 0,
                   s_id = 0,
                   txt = paste0("\n\n## ", unique(.$text_title), "\n\n"),
                   cf = NA,
                   hn = NA,
                   msa = NA,
                   gls = NA,
                   free_trans = NA), .before = 1)  %>% 
    arrange(text_title, p_id, s_id) %>% 
    select(-c(all_rows[!(all_rows %in% params$rows)])) %>% 
    pivot_longer(names_to = "type", values_to = "value", cols = c("txt", params$rows, "free_trans")) %>% 
    na.omit() %>% 
    ungroup() %>% 
    select(value) %>% 
    unlist() %>% 
    unname() %>% 
    cat(sep = "")  
}
```
---
title: "Phonetic fieldwork and experiments with the `phonfieldwork` package for R"
author: "George Moroz"
institute: "Linguistic Convergence Laboratory, NRU HSE"
date: |
    | 28 September 2020
    |
    | Presentation is available here: \alert{tinyurl.com/y2x8ppnl}
    | ![](images/00_qrcode.png)'
output: 
  beamer_presentation:
    latex_engine: xelatex
    citation_package: natbib
    keep_tex: false
    includes:
      in_header: "config/presento.sty"
bibliography: bibliography.bib
biblio-style: "apalike"
editor_options: 
  chunk_output_type: console
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE, warning = FALSE, message = FALSE)
# library(qrcode)
# png(filename="images/00_qrcode.png", width = 200, height = 200)
# qrcode_gen("https://github.com/agricolamz/phonfieldwork/raw/master/presentations/2020.29.09_Conlab_phonfieldwork_talk/2020.29.09_Conlab_phonfieldwork_talk.pdf")
# dev.off()

unlink("data", recursive = TRUE)
dir.create("data")
dir.create("data/s1")
dir.create("data/s2")
file.copy(paste0("../backup/s1/", list.files("../backup/s1/")), "data/s1/")
file.copy(paste0("../backup/s2/", list.files("../backup/s2/")), "data/s2/")
file.copy(paste0("../backup/", list.files("../backup/", pattern = "test")), "data/")
```

# About `phonfieldwork`

## About `phonfieldwork`

This package was started as a help for our guest student Margaux Dubuis from Switzerland, who was going to study New Caledonian language Vamale (`vama1243`).

```{r, echo=FALSE, out.width = '60%'}
knitr::include_graphics("images/01_Margaux.jpg")
```

## About `phonfieldwork`

* 2019-08-24 --- first release on CRAN
* Now there is a v0.0.7 on CRAN, and total 6363 downloads

```{r, echo=FALSE, out.width = '100%'}
knitr::include_graphics("images/02_commits.png")
```

## About `phonfieldwork`

* On June 20 I have made [a submission](https://github.com/ropensci/software-review/issues/385) to [rOpenSci](https://ropensci.org/).
* Here are my reviewers:

```{r, echo=FALSE, out.width = '90%'}
knitr::include_graphics("images/03_Niko_Jonathan.png")
```

# Introduction

## Most phonetic research consists of the following steps:

1. Formulate a research question. Think of what kind of data is necessary to answer this question, what is the appropriate amount of data, what kind of annotation you will do, what kind of statistical models and visualizations you will use, etc.
2. Create a list of stimuli.
3. Elicite list of stimuli from speakers who signed an Informed Consent statement, agreeing to participate in the experiment to be recorded on audio and/or video.
4. Annotate the collected data.
5. Extract the information from annotated data.
6. Create visualizations.
7. Evaluate your statistical models.
8. Report your results.
9. Publish your data. \pause

The `phonfieldwork` package is created for helping with items 3, 4 (partially), 5, 6, and 9.

## Why/when do you need the phonfieldwork package?

These ideal plan hides some technical subtasks:

* creating a presentation for elicitation task
* renaming and concatenating multiple sound files
* automatic annotation in Praat TextGrids [@boersma19]
* creating a searchable `.html` table with annotations, spectrograms and ability to hear sound
* converting multiple formats (Praat, ELAN [@brugman04] and EXMARaLDA [@schmidt09], .srt and some others) \pause

All of these tasks can be solved by a mixture of different tools:

* any programming language can handle automatic file renaming
* Praat contains scripts for concatenating and renaming files \pause

`phonfieldwork` provides a functionality that will make it easier to solve those tasks independently of any additional tools. You can also compare the functionality with other packages: `rPraat` [@borzil16], `textgRid` [@reidy16], `pympi` [@lubbers13]

## Philosophy of the `phonfieldwork` package

* each stimulus is a separate file
* researcher carefully listens to consultants to make sure that they are producing the kind of speech they wanted
* in case a speaker does not produce clear repetitions, researcher ask them to repeat the task

There are some phoneticians who prefer to record everything for language documentation purposes. I think that should be a separate task. If you insist on recording everything, it is possible to run two recorders at the same time: one could run during the whole session, while the other is used to produce small audio files. You can also use special software to record your stimuli automatically on a computer (e.g. [`SpeechRecorder`](https://www.bas.uni-muenchen.de/Bas/software/speechrecorder/) [@draxler04],  `PsychoPy` [@peirce19]).
   
# Installation of the package

## Install `phonfieldwork`

`phonfieldwork` is an R package, so you need to install [R](https://cloud.r-project.org/), [RStudio](https://rstudio.com/products/rstudio/download/#download) (optional) or use [rstudio.cloud](rstudio.cloud). There are two possibilities for installing package in R:

* official version from CRAN

```{r, eval = FALSE}
install.packages("phonfieldwork")
```

* development version from GitHub

```{r, eval = FALSE}
devtools::install_github("agricolamz/phonfieldwork")
```

\pause

Since this package is under [rOpenSci review](https://github.com/ropensci/software-review/issues/385) there is a chance that in couple months the adress could be changed to `ropensci/phonfieldwork`, and documentation page will be moved from [agricolamz.github.io/phonfieldwork](agricolamz.github.io/phonfieldwork) to
[docs.ropensci.org/phonfieldwork](docs.ropensci.org/phonfieldwork).

```{r}
library("phonfieldwork")
packageVersion("phonfieldwork") # Unreleased version
```

# Creating your presentation

##   Make a list of your stimuli

There are several ways to enter information about a list of stimuli into R:

* list them with the `c()` command

```{r}
my_stimuli <- c("tip", "tap", "top")
```

* read from `.csv` file with the `read.csv()` function:

```{r, eval = FALSE}
my_stimuli_df <- read.csv("my_stimuli_df.csv")
```

* read from `.xls` or `.xlsx` file with the `read_xls()` or `read_xlsx()` function from the `readxl` package (run `install.packages("readxl")` in case you don't have it installed):

```{r, eval = FALSE}
library("readxl")
my_stimuli_df <- read_xlsx("my_stimuli_df.xlsx")
```


##  Create a presentation based on a list of stimuli

```{r}
create_presentation(stimuli = my_stimuli, # it is "tip" "tap" "top"
                    output_file = "first_example",
                    output_dir = "data/")
```

Here is [the result](https://agricolamz.github.io/phonfieldwork/additional/first_example.html). \pause

It is also possible to use images (or `.gif`) as a stimuli:

```{r}
my_image <- system.file("extdata", "r-logo.png", package = "phonfieldwork")
my_image

create_presentation(stimuli = c("rzeka", "drzewo", my_image),
                    external = 3,
                    output_file = "second_example",
                    output_dir = "data/")
```

Here is [the result](https://agricolamz.github.io/phonfieldwork/additional/second_example.html).

# Data renaming

## Obtained data

After collecting data and removing soundfiles with unsuccesful elicitations, one could end up with the following structure:

```{bash, echo = FALSE}
tree data | tail -n 13 | head -n 8
```

## Rename collected data

```{r}
rename_soundfiles(stimuli = my_stimuli, # it is "tip" "tap" "top"
                  prefix = "s1_",
                  path = "data/s1/")
```

```{bash, echo = FALSE}
tree data | tail -n 18 | head -n 13
```

## Rename collected data

Here is the contents of `logging.csv`:

```{r}
read.csv("data/s1/backup/logging.csv")
```

Niko Partanen also nicely suggested to have a `missing` argument for `rename_soundfiles()` in order to make handling missing soundfiles easier.

## Rename collected data

```{r}
rename_soundfiles(stimuli = my_stimuli, # it is "tip" "tap" "top"
                  prefix = paste0(1:3, "_"),
                  suffix = "_s2",
                  path = "data/s2/",
                  backup = FALSE,
                  logging = FALSE,
                  autonumbering = FALSE)
```

```{bash, echo = FALSE}
tree data | tail -n 9 | head -n 4
```

## Get sound duration

Sometimes it is useful to get information about sound duration:

```{r}
get_sound_duration("data/s1/2_s1_tap.wav")
```

```{r}
get_sound_duration(sounds_from_folder = "data/s2/")
```

# Data merging

## Merge all data together

After all the files are renamed, you can merge them into one:

```{r}
concatenate_soundfiles(path = "data/s1/",
                       result_file_name = "s1_all")
```

```{bash, echo = FALSE}
tree data/s1 | head -n 11
```

## Merge all data together

```{r, echo=FALSE}
draw_sound("data/s1/s1_all.wav", "data/s1/s1_all.TextGrid", title = "s1_all")
```

# Data annotation

## Annotate your data

It is possible to create annotation in advance:

```{r}
annotate_textgrid(annotation = my_stimuli,
                  textgrid = "data/s1/s1_all.TextGrid")
```

## Annotate your data

```{r, echo=FALSE}
draw_sound("data/s1/s1_all.wav", "data/s1/s1_all.TextGrid", title = "s1_all")
```

## Create subannotation

Imagine that we are interested in annotation of vowels. The most common solution will be open Praat and create new annotations. But it is also possible to create them in advance. The idea that you choose some baseline tier that later will be automatically cutted into smaller pieces on the other tier.

```{r}
create_subannotation(textgrid = "data/s1/s1_all.TextGrid",
                     tier = 1, # this is a baseline tier
                     n_of_annotations = 3) # how many empty annotations per unit?
```

## Annotate subannotation in advance

```{r, echo=FALSE}
draw_sound("data/s1/s1_all.wav", "data/s1/s1_all.TextGrid", title = "s1_all")
```

## Annotate subannotation in advance

Now we can annotate created tier:
```{r}
annotate_textgrid(annotation = c("", "ı", "", "", "æ", "", "", "ɒ", ""),
                  textgrid = "data/s1/s1_all.TextGrid",
                  tier = 3,
                  backup = FALSE)
```

List annotations by hand is a boring task, so if you have a prepared list of annotations, the merege could be done with the following code:
```{r}
vowels <- c("ı", "æ", "ɒ")
unlist(lapply(vowels, function(x){c("", x, "")}))
```


## Create subannotation

```{r, echo=FALSE, eval=FALSE}
draw_sound("data/s1/s1_all.wav", "data/s1/s1_all.TextGrid", title = "s1_all", output_file = "images/04_subannotation")
```

```{r, echo=FALSE, out.width = '100%'}
knitr::include_graphics("images/04_subannotation.png")
```

## The only thing left is to move annotation boundaries

```{r, include=FALSE}
file.copy("../backup/s1_all.TextGrid", to = "data/s1/s1_all.TextGrid", overwrite = TRUE)
```

```{r, echo=FALSE, eval=FALSE}
draw_sound("data/s1/s1_all.wav", "data/s1/s1_all.TextGrid", title = "s1_all", output_file = "images/05_subannotation")
```

```{r, echo=FALSE, out.width = '100%'}
knitr::include_graphics("images/05_subannotation.png")
```

# Data extraction
## Data viewer

Sound viewer (here is an [example](https://agricolamz.github.io/phonfieldwork/additional/stimuli_viewer.html)) is a useful tool that combine together your annotations and make it searchable. It is also produce a ready to go .html file that could be uploaded on the server (e. g. to Github Pages) and be availible for anyone in the world.

In order to do that we need:

* seperate folder with soundfiles
* separate folder with spectorgrams (optional)
* `data.frame` with data about utterances or speakers

## Data extraction
First, it is important to create a folder where all of the extracted files will be stored:

```{r}
dir.create("data/s1/s1_sounds")
```

It is possible extract to extract all annotated files based on an annotation tier:

```{r}
extract_intervals(file_name = "data/s1/s1_all.wav",
                  textgrid = "data/s1/s1_all.TextGrid",
                  tier = 3,
                  path = "data/s1/s1_sounds/",
                  prefix = "s1_")
```

## Data extraction

After those commands, one could end up with the following structure:

```{bash, echo = FALSE}
tree data/s1 | head -n 15
```

# Data visulization

## Sound visulization in `phonfieldwork`

The easiest way to visualise sound in phonfieldwork:
```{r, fig.height=6}
draw_sound(file_name = "data/s1/2_s1_tap.wav")
```

## Sound visulization in `phonfieldwork`

```{r, eval=FALSE}
draw_sound(file_name = "data/s1/s1_all.wav", 
           annotation = "data/s1/s1_all.TextGrid")
```

```{r, include=FALSE, eval=FALSE}
draw_sound(file_name = "data/s1/s1_all.wav", 
           annotation = "data/s1/s1_all.TextGrid", 
           output_file = "images/06_sound_visualisation")
```

```{r, echo=FALSE, out.width = '100%'}
knitr::include_graphics("images/06_sound_visualisation.png")
```

```{r, eval = FALSE}
draw_sound("data/s1/s1_all.wav",
           "data/s1/s1_all.TextGrid",
           from = 0.4, to = 0.95)
```

```{r, include=FALSE, eval = FALSE}
draw_sound(file_name = "data/s1/s1_all.wav", 
           annotation = "data/s1/s1_all.TextGrid", 
           output_file = "images/07_sound_visualisation")
```


```{r, echo=FALSE, out.width = '100%'}
knitr::include_graphics("images/07_sound_visualisation.png")
```


## Create multiple spectrograms

```{r}
draw_sound(sounds_from_folder = "data/s1/s1_sounds/",
           pic_folder_name = "s1_pics")
```

```{bash, echo = FALSE}
tree data/s1 head -n 18
```

# Creating a data viewer -- template for the data sharing

## Create a [viewer](https://agricolamz.github.io/phonfieldwork/additional/stimuli_viewer.html)

* seperate folder with soundfiles 
* separate folder with spectorgrams (optional) 
* \alert{→ data.frame with data about utterances or speakers}

```{r}
df <- data.frame(word  = c("tap", "tip", "top"),
                 sounds = c("æ", "ı", "ɒ"))
df

create_viewer(audio_dir = "data/s1/s1_sounds/",
              picture_dir = "data/s1/s1_pics/",
              table = df,
              output_dir = "data/s1/",
              output_file = "stimuli_viewer")
```

# Reading data from different linguistic sources

##  Read linguistic files into R

* `textgrid_to_df()` (Praat)
* `eaf_to_df()` (ELAN)
* `exb_to_df()` (EXMARaLDA)
* `srt_to_df()` (subtitle file)
* `audacity_to_df()` (Audacity)
* `flextext_to_df()` (FieldWorks)

##  Read linguistic files into R

```{r}
draw_sound(file_name = "data/test.wav",
           annotation = eaf_to_df("data/test.eaf"))
```


# Get help and cite

## Get help and cite

You can always write an email or open an [issue on GitHub](https://github.com/agricolamz/phonfieldwork/issues), asking some questions.

The most recent citation information is avalible with this command:

```{r}
citation("phonfieldwork")
```

# References {.allowframebreaks}
---
title: "Phonetic Fieldwork and Experiments with the phonfieldwork Package for R"
author: 
  - George Moroz
date: "2020-12-22"
slug: phonfieldwork-phonetic-fieldwork-and-experiments
categories:
  - blog
tags:
  - Software Peer Review
  - packages
  - R
  - linguistics
  - phonetics
  - sound analysis
  - fieldwork
  - phonfieldwork
  - lingtypology
  - community
package_version: 0.0.10
description: "phonfieldwork package provides a huge variety of tools that can be used by phonetic researchers in the field and in laboratories."
output:
  html_document:
    keep_md: yes
---

```{r setup, include=FALSE}
# Options to have images saved in the post folder
# And to disable symbols before output
knitr::opts_chunk$set(fig.path = "", comment = "")
# knitr hook to make images output use Hugo options
knitr::knit_hooks$set(
  plot = function(x, options) {
    hugoopts <- options$hugoopts
    paste0(
      "{{<figure src=",
      '"', x, '" ',
      if (!is.null(hugoopts)) {
        glue::glue_collapse(
          glue::glue('{names(hugoopts)}="{hugoopts}"'),
          sep = " "
        )
      },
      ">}}\n"
    )
  }
)
# knitr hook to use Hugo highlighting options
knitr::knit_hooks$set(
  source = function(x, options) {
    hlopts <- options$hlopts
    paste0(
      "```r ",
      if (!is.null(hlopts)) {
        paste0("{",
               glue::glue_collapse(
                 glue::glue('{names(hlopts)}={hlopts}'),
                 sep = ","
               ), "}"
        )
      },
      "\n", glue::glue_collapse(x, sep = "\n"), "\n```\n"
    )
  }
)
```

## Science craft

As a field linguist, I have spent a lot of time working in villages in the Caucasus, collecting audio from speakers of indigenous languages. The processing of such data involves a lot of time-consuming tasks, so during my field trips I created my own pipeline for data collection. I wrote a number of scripts using different programming languages for automatic renaming, merging, and preannotation of files, making backups, visualizing some data, etc. My method consisted of a combination of solutions developed on the fly to solve specific, independent tasks, without thinking of the future.

This became a problem when I wanted to pass on my knowledge to my students. Not all of them were familiar with all of the programming languages I used, some of the code had become outdated, and some code did not work properly on all operating systems.

The practical side of linguistic research is often neglected. Linguists rely on manual and makeshift solutions for working with data and simultaneously write similar scripts, reinventing the same wheel as part of the process. The few existing software programs are undercited: there are a lot of papers presenting data processed in R,[^1] Praat,[^2] ELAN[^3] or other software, which do not cite them. Methodologically oriented papers are valued mainly for their theoretical implications.

Recent years have seen increased attention to best practices in linguistic data management, for example through the introduction of [Cross-Linguistic Data Formats](https://cldf.clld.org/).[^4] I think the time is also ripe to pay more attention to what I like to call **science craft** -- improving the methods for curating linguistic data, and allowing future generations of linguists to benefit from our experience.

That is why I decided to create a toolkit for phonetic researchers in the form of an R package, [phonfieldwork](https://docs.ropensci.org/phonfieldwork): it is written using one programming language, and easy to use for non-coders and people who are not familiar with R.

Below is an overview of its main functions.

## phonfieldwork

Most phonetic research consists of the following steps:

1. Formulate a research question. Think of what kind of data is necessary to answer this question, what is the appropriate amount of data, what kind of annotation you will do, what kind of statistical models and visualizations you will use, etc.
2. Create a list of stimuli.
3. Elicit the list of stimuli from speakers who signed an Informed Consent statement, agreeing to participate in the experiment and to be recorded on audio and/or video. Keep an eye on the recording settings: the sampling rate, resolution (bit), and number of channels should be the same across all recordings.
4. Annotate the collected data.
5. Extract the collected data.
6. Create visualizations and evaluate your statistical models.
7. Report your results.
8. Publish your data.

The phonfieldwork package is created for helping with items 3, 5, and 8, and partially 4. If you are interested in the whole pipeline please read the whole [Get started section](https://docs.ropensci.org/phonfieldwork/articles/phonfieldwork.html) in the package documentation. Below I will introduce some specific features.

##  What can be done with phonfieldwork?

First let's load the package:

```{r load_the_library}
library(phonfieldwork)
```

### Create a presentation based on a list of stimuli

It is easier to collect your stimuli with a presentation. In order to create it you need a list of stimuli (as an example I will use a simple vector, but of course it can be a column from `.csv` or `.xlsx` files):

```{r, eval=FALSE}
create_presentation(stimuli = c("tip", "tap", "top"),
                    output_file = "first_example",
                    output_dir = "...path/to/your/folder")
```

As a result, a file `first_example.html` will appear in the output folder ([here](https://ropensci.github.io/phonfieldwork/additional/first_example.html) is an example of a such file[^7]). It is also possible to add translations and even use images (see [the documentation page](https://docs.ropensci.org/phonfieldwork/articles/phonfieldwork.html#create-a-presentation-based-on-a-list-of-stimuli-1)).

### Sound annotation formats

There was a goal to be able to convert multiple sound annotation formats into data.frame format, so you can find a whole bunch of functions that serve this purpose:

* convert Praat `.TextGrid` files[^6]; see also [`rPraat`](https://fu.ff.cuni.cz/praat/#rpraat-package-for-r) and [`textgRid`](https://github.com/patrickreidy/textgRid) packages
```{r textgrid_to_df_example}
textgrid_to_df(system.file("extdata", "test.TextGrid", package = "phonfieldwork"))
```

* convert ELAN `.eaf` files; see also the [FRelan](https://github.com/langdoc/FRelan) package by Niko Partanen

```{r eaf_to_df_example}
eaf_to_df(system.file("extdata", "test.eaf", package = "phonfieldwork"))
```

* convert EXMARaLDA `.exb` files

```{r exb_to_df_example}
exb_to_df(system.file("extdata", "test.exb", package = "phonfieldwork"))
```

* convert subtitles `.srt` file

```{r srt_to_df_example}
srt_to_df(system.file("extdata", "test.srt", package = "phonfieldwork"))
```

* convert Audacity `.txt` file

```{r audacity_to_df_example}
audacity_to_df(system.file("extdata", "test_audacity.txt", package = "phonfieldwork"))
```

There is also an option to work with `.flextext` files from FLEx, but since this is not relevant for phonetics, I will skip this part.

### Sound preannotation

The sound annotation process takes a lot of time and in most cases it is really boring. Because of this boredom annotators make mistakes during the annotation process. In order to prevent it phonfieldwork provides a bunch of functions that preannotate your sound files (`annotate_textgrid()`, `create_subannotation()` and  `create_empty_textgrid()`). A detailed example would be too large to put in a blog post (you can find it in the [documentation](https://docs.ropensci.org/phonfieldwork/articles/phonfieldwork.html#annotate-your-data)). Here is a `.gif` (created with [magick](https://docs.ropensci.org/magick/index.html)) that illustrates the whole process:

```{r create_gif, echo=FALSE, eval=FALSE}
library(magick)
im1 <- image_read("01.png")
im2 <- image_read("02.png")
im3 <- image_read("03.png")
im4 <- image_read("04.png")
image_animate(image_scale(c(im1, im2, im3, im4), "700"), fps = 0.5, dispose = "previous")
```

<!--html_preserve-->
{{< figure src = "annotation_in_phonfieldwork.gif" width = "700" alt = "Oscillogram and spectrogram with text annotations below. GIF shows how phonfieldwork can add preannotations to a recording" >}}
<!--/html_preserve-->

As you can see, it is possible to create all the annotations in advance and leave the annotators with one task only: to move boundaries.

### Sound visualization

Sound visualization is a common task that could be solved via different programs and R packages, but it can also be done with phonfieldwork:

```{r spec_example, hugoopts=list(alt="Example of an oscillogram and a spectrogram", caption="Example of an oscillogram and a spectrogram generated by the draw_sound() function", width=600)}
file <- system.file("extdata", "test.wav", package = "phonfieldwork")
draw_sound(file)
```

A feature specific to phonfieldwork, is that it is possible to zoom in to some part of the sound with a spectrogram (time in the `zoom` argument is represented in seconds):

```{r spec_zoom_example, hugoopts=list(alt="Oscillogram and spectrogram zoomed in to a specific time to illustrate usage of the zoom argument", caption="Usage of the zoom argument", width=600)}
draw_sound(file, zoom = c(0.2, 0.4))
```

It is also possible to visualize any sound annotation format that was converted to dataframe:

```{r spec_annotation_example, hugoopts=list(alt="Oscillogram and spectrogram with text annotations below the figures, illustrating usage of the annotation argument", caption="Usage of the annotation argument", width=600)}
our_textgrid <- system.file("extdata", "test.TextGrid", package = "phonfieldwork")
draw_sound(file, 
           annotation = textgrid_to_df(our_textgrid))
```

### Sound viewer

If you have folders with small sound chunks and their visualizations it is possible to create a sound viewer like [this](https://ropensci.github.io/phonfieldwork/additional/stimuli_viewer.html). This is done using the `create_viewer()` function:

```{r viewer_example, eval = FALSE}
create_viewer(audio_dir = ".../sounds/", # path to folder with sounds
              picture_dir = ".../pictures/", # path to folder with pictures
              table = df, # dataframe with additional information
              output_dir = "...", # where to store the result?
              output_file = "...") # how to name the result file?
```

If you are familiar with my package lingtypology[^5] for interactive linguistic map generation and API for typological databases, there is good news for you: it is possible to connect the two packages, creating an interactive map that shares the same hear and view buttons. Here is [an example](https://ropensci.github.io/phonfieldwork/additional/stimuli_viewer2.html).

I really hope that this format will become a new tool for searching, analyzing and sharing phonetic data. However, there is always a risk that this tool can be misused, so please read the [text about Ethical Research with phonfieldwork](https://docs.ropensci.org/phonfieldwork/articles/ethical_research_with_phonfieldwork.html).

It is only a brief introduction into the phonfieldwork functionality, for more details read the whole [Get started section](https://docs.ropensci.org/phonfieldwork/articles/phonfieldwork.html) in the package documentation.

## Acknowledgements

I would like to thank

* [my rOpenSci reviewers](https://github.com/ropensci/software-review/issues/385) [Jonathan Keane](/author/jonathan-keane/) and [Niko Partanen](https://github.com/nikopartanen) for their interesting comments and ideas;
* [Melina Vidoni](/author/melina-vidoni/) for being the package review editor;
* participants of seminars at the School of linguistics and Linguistic Convergence Laboratory at HSE, Moscow, where I first presented this package;
* my friends [Vanya Kapitonov](https://github.com/vkusvody), [Neige Rochant](https://github.com/Neigelily), and [Samira Verhees](https://github.com/sverhees/) for sharing their data and problems with me, which has made phonfieldwork better.

[^1]: R Core Team. 2020. R: A Language and Environment for Statistical Computing. Vienna, Austria: R Foundation for Statistical Computing. <https://www.R-project.org/>.
[^2]: Boersma, P., and D. Weenink. 2020. "Praat: Doing Phonetics by Computer [Computer Program]. Version 6.1.35, Retrieved 5 December 2020 from https://Www.praat.org/."
[^3]: Brugman, Hennie, Albert Russel, and Xd Nijmegen. 2004. "Annotating Multi-Media/Multi-Modal Resources with ELAN." In LREC.
[^4]: Forkel, Robert, Johann-Mattis List, Simon J. Greenhill, Christoph Rzymski, Sebastian Bank, Harald Hammarström, Martin Haspelmath, and Russell D. Gray. 2018. “Cross-Linguistic Data Formats, Advancing Data Sharing and Re-Use in Comparative Linguistics.” Sci Data 5 (1): 1–10.
[^5]: Moroz, George. 2017. Lingtypology: Easy Mapping for Linguistic Typology. https://CRAN.R-project.org/package=lingtypology.
[^6]: Just change the `system.file()` function to the path to the file.
[^7]: It will not look nice on a phone, but you can use it on a PC or a tablet.---
title: "Phonetic fieldwork and experiments with the `phonfieldwork` package for R"
author: "George Moroz"
institute: "Linguistic Convergence Laboratory, NRU HSE"
date: |
    | 08 August 2020, Grupo de estadística para el estudio del lenguaje
    |
    |
    | Presentation is available here: \alert{tinyurl.com/y6lf5ch4}
    | ![](images/01_qrcode.png)'
output: 
  beamer_presentation:
    latex_engine: xelatex
    citation_package: natbib
    keep_tex: false
    includes:
      in_header: "config/presento.sty"
bibliography: bibliography.bib
biblio-style: "apalike"
editor_options: 
  chunk_output_type: console
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE, warning = FALSE, message = FALSE)
# library(qrcode)
# png(filename="images/01_qrcode.png", width = 200, height = 200)
# qrcode_gen("https://github.com/agricolamz/phonfieldwork/raw/master/presentations/2020.08.08_Grupo_de_estad%C3%ADstica_para_el_estudio_del_lenguaje/2020.08.08_Grupo_de_estad%C3%ADstica_para_el_estudio_del_lenguaje_moroz_phonfieldwork.pdf")
# dev.off()

unlink("data", recursive = TRUE)
dir.create("data")
dir.create("data/s1")
dir.create("data/s2")
file.copy(paste0("../backup/s1/", list.files("../backup/s1/")), "data/s1/")
file.copy(paste0("../backup/s2/", list.files("../backup/s2/")), "data/s2/")
file.copy(paste0("../backup/", list.files("../backup/", pattern = "test")), "data/")
```

# Introduction

## Most phonetic research consists of the following steps:

1. Formulate a research question. Think of what kind of data is necessary to answer this question, what is the appropriate amount of data, what kind of annotation you will do, what kind of statistical models and visualizations you will use, etc.
2. Create a list of stimuli.
3. Elicite list of stimuli from speakers who signed an Informed Consent statement, agreeing to participate in the experiment to be recorded on audio and/or video.
4. Annotate the collected data.
5. Extract the information from annotated data.
6. Create visualizations.
7. Evaluate your statistical models.
8. Report your results.
9. Publish your data. \pause

The `phonfieldwork` package is created for helping with items 3, 4 (partially), 5, 6, and 9.

## Why/when do you need the phonfieldwork package?

These ideal plan hides some technical subtasks:

* creating a presentation for elicitation task
* renaming and concatenating multiple sound files
* automatic annotation in Praat TextGrids [@boersma19]
* creating a searchable `.html` table with annotations, spectrograms and ability to hear sound
* converting multiple formats (Praat, ELAN [@brugman04] and EXMARaLDA [@schmidt09]) \pause

All of these tasks can be solved by a mixture of different tools:

* any programming language can handle automatic file renaming
* Praat contains scripts for concatenating and renaming files \pause

`phonfieldwork` provides a functionality that will make it easier to solve those tasks independently of any additional tools. You can also compare the functionality with other packages: `rPraat` [@borzil16], `textgRid` [@reidy16], `pympi` [@lubbers13]

## Philosophy of the `phonfieldwork` package

* each stimulus is a separate file
* researcher carefully listens to consultants to make sure that they are producing the kind of speech they wanted
* in case a speaker does not produce clear repetitions, researcher ask them to repeat the task

There are some phoneticians who prefer to record everything for language documentation purposes. I think that should be a separate task. If you insist on recording everything, it is possible to run two recorders at the same time: one could run during the whole session, while the other is used to produce small audio files. You can also use special software to record your stimuli automatically on a computer (e.g. `PsychoPy` [@peirce19]).
   
# Installation of the package

## Install `phonfieldwork`

`phonfieldwork` is an R package, so you need to install [R](https://cloud.r-project.org/), [RStudio](https://rstudio.com/products/rstudio/download/#download) (optional) or use [rstudio.cloud](rstudio.cloud). There are two possibilities for installing package in R:

* official version from CRAN

```{r, eval = FALSE}
install.packages("phonfieldwork")
```

* development version from GitHub

```{r, eval = FALSE}
devtools::install_github("agricolamz/phonfieldwork")
```

\pause

Since this package is under [rOpenSci review](https://github.com/ropensci/software-review/issues/385) there is a chance that in couple months the adress could be changed to `ropensci/phonfieldwork`, and documentation page will be moved from [agricolamz.github.io/phonfieldwork](agricolamz.github.io/phonfieldwork) to
[docs.ropensci.org/phonfieldwork](docs.ropensci.org/phonfieldwork).

```{r}
library("phonfieldwork")
packageVersion("phonfieldwork") # Unreleased version
```

# Creating your presentation

##   Make a list of your stimuli

There are several ways to enter information about a list of stimuli into R:

* list them with the `c()` command

```{r}
my_stimuli <- c("tip", "tap", "top")
```

* read from `.csv` file with the `read.csv()` function:

```{r, eval = FALSE}
my_stimuli_df <- read.csv("my_stimuli_df.csv")
```

* read from `.xls` or `.xlsx` file with the `read_xls()` or `read_xlsx()` function from the `readxl` package (run `install.packages("readxl")` in case you don't have it installed):

```{r, eval = FALSE}
library("readxl")
my_stimuli_df <- read_xlsx("my_stimuli_df.xlsx")
```


##  Create a presentation based on a list of stimuli

```{r}
create_presentation(stimuli = my_stimuli, # it is "tip" "tap" "top"
                    output_file = "first_example",
                    output_dir = "data/")
```

Here is [the result](https://agricolamz.github.io/phonfieldwork/additional/first_example.html). \pause

It is also possible to use images (or `.gif`) as a stimuli:

```{r}
my_image <- system.file("extdata", "r-logo.png", package = "phonfieldwork")
my_image

create_presentation(stimuli = c("rzeka", "drzewo", my_image),
                    external = 3,
                    output_file = "second_example",
                    output_dir = "data/")
```

Here is [the result](https://agricolamz.github.io/phonfieldwork/additional/second_example.html).

# Data renaming

## Obtained data

After collecting data and removing soundfiles with unsuccesful elicitations, one could end up with the following structure:

```{bash, echo = FALSE}
tree data | tail -n 13 | head -n 8
```

## Rename collected data

```{r}
rename_soundfiles(stimuli = my_stimuli, # it is "tip" "tap" "top"
                  prefix = "s1_",
                  path = "data/s1/")
```

```{bash, echo = FALSE}
tree data | tail -n 17 | head -n 12
```

## Rename collected data

```{r}
rename_soundfiles(stimuli = my_stimuli, # it is "tip" "tap" "top"
                  prefix = paste0(1:3, "_"),
                  suffix = "_s2",
                  path = "data/s2/",
                  backup = FALSE)
```

```{bash, echo = FALSE}
tree data | tail -n 17 | head -n 12
```

## Rename collected data

Sometimes it is better to keep and order that will deal with the computer sorting:

```{r}
add_leading_symbols(1:105)
```

So it is better to use the result of `add_leading_symbols()` as a prefix during the renaming of a huge amount of files.

## Get sound duration

Sometimes it is useful to get information about sound duration:

```{r}
get_sound_duration("data/s1/s1_tap.wav")
```

```{r}
get_sound_duration(sounds_from_folder = "data/s2/")
```

# Data merging

## Merge all data together

After all the files are renamed, you can merge them into one:

```{r}
concatenate_soundfiles(path = "data/s1/",
                       result_file_name = "s1_all")
```

```{bash, echo = FALSE}
tree data/s1 | head -n 10
```

## Merge all data together

```{r, echo=FALSE}
draw_sound("data/s1/s1_all.wav", "data/s1/s1_all.TextGrid", title = "s1_all")
```

# Data annotation

## Annotate your data

It is possible to create annotation in advance (since file concatination is made according to files sorted on the comuter I use the `sort()` function in order to make correct annotation):

```{r}
annotate_textgrid(annotation =  sort(my_stimuli),
                  textgrid = "data/s1/s1_all.TextGrid")
```

## Annotate your data

```{r, echo=FALSE}
draw_sound("data/s1/s1_all.wav", "data/s1/s1_all.TextGrid", title = "s1_all")
```

## Create subannotation

Imagine that we are interested in annotation of vowels. The most common solution will be open Praat and create new annotations. But it is also possible to create them in advance. The idea that you choose some baseline tier that later will be automatically cutted into smaller pieces on the other tier.

```{r}
create_subannotation(textgrid = "data/s1/s1_all.TextGrid",
                     tier = 1, # this is a baseline tier
                     n_of_annotations = 3) # how many empty annotations per unit?
```

## Annotate subannotation in advance

```{r, echo=FALSE}
draw_sound("data/s1/s1_all.wav", "data/s1/s1_all.TextGrid", title = "s1_all")
```

## Annotate subannotation in advance

Now we can annotate created tier:
```{r}
annotate_textgrid(annotation = c("", "æ", "", "", "ı", "", "", "ɒ", ""),
                  textgrid = "data/s1/s1_all.TextGrid",
                  tier = 3,
                  backup = FALSE)
```

List annotations by hand is a boring task, so if you have a prepared list of annotations, the merege could be done with the following code:
```{r}
vowels <- c("æ", "ı", "ɒ")
unlist(lapply(vowels, function(x){c("", x, "")}))
```


## Create subannotation

```{r, echo=FALSE, eval = FALSE}
draw_sound("data/s1/s1_all.wav", "data/s1/s1_all.TextGrid", title = "s1_all", output_file = "images/02_subannotation")
```

![](images/02_subannotation.png)

## The only thing left is to move annotation boundaries

```{r, include=FALSE}
file.copy("../backup/s1_all.TextGrid", to = "data/s1/s1_all.TextGrid", overwrite = TRUE)
```

```{r, echo=FALSE, eval=FALSE}
draw_sound("data/s1/s1_all.wav", "data/s1/s1_all.TextGrid", title = "s1_all", output_file = "images/03_subannotation")
```

![](images/03_subannotation.png)

# Data extraction
## Data viewer

Sound viewer (here is an [example](https://agricolamz.github.io/phonfieldwork/additional/stimuli_viewer.html)) is a useful tool that combine together your annotations and make it searchable. It is also produce a ready to go .html file that could be uploaded on the server (e. g. to Github Pages) and be availible for anyone in the world.

In order to do that we need:

* seperate folder with soundfiles
* separate folder with spectorgrams (optional)
* `data.frame` with data about utterances or speakers

## Data extraction
First, it is important to create a folder where all of the extracted files will be stored:

```{r}
dir.create("data/s1/s1_sounds")
```

It is possible extract to extract all annotated files based on an annotation tier:

```{r}
extract_intervals(file_name = "data/s1/s1_all.wav",
                  textgrid = "data/s1/s1_all.TextGrid",
                  tier = 3,
                  path = "data/s1/s1_sounds/",
                  prefix = "s1_")
```

## Data extraction

After those commands, one could end up with the following structure:

```{bash, echo = FALSE}
tree data/s1 | head -n 14
```

# Data visulization

## Sound visulization in `phonfieldwork`

The easiest way to visualise sound in phonfieldwork:
```{r}
draw_sound(file_name = "data/s1/s1_tap.wav")
```

## Sound visulization in `phonfieldwork`

```{r, eval=FALSE}
draw_sound(file_name = "data/s1/s1_all.wav", 
           annotation = "data/s1/s1_all.TextGrid")
```

```{r, include=FALSE, eval=FALSE}
draw_sound(file_name = "data/s1/s1_all.wav", 
           annotation = "data/s1/s1_all.TextGrid", 
           output_file = "images/04_sound_visualisation.png")
```

![](images/04_sound_visualisation.png)


```{r, eval = FALSE}
draw_sound("data/s1/s1_all.wav",
           "data/s1/s1_all.TextGrid",
           from = 0.4, to = 0.95)
```

```{r, include=FALSE, eval = FALSE}
draw_sound(file_name = "data/s1/s1_all.wav", 
           annotation = "data/s1/s1_all.TextGrid", 
           output_file = "images/05_sound_visualisation.png")
```

![](images/05_sound_visualisation.png)

## Create multiple spectrograms

```{r}
draw_sound(sounds_from_folder = "data/s1/s1_sounds/",
           pic_folder_name = "s1_pics")
```

```{bash, echo = FALSE}
tree data/s1 head -n 18
```

# Creating a data viewer -- template for the data sharing

## Create a [viewer](https://agricolamz.github.io/phonfieldwork/additional/stimuli_viewer.html)

* seperate folder with soundfiles 
* separate folder with spectorgrams (optional) 
* \alert{→ data.frame with data about utterances or speakers}

```{r}
df <- data.frame(word  = c("tap", "tip", "top"),
                 sounds = c("æ", "ı", "ɒ"))
df

create_viewer(audio_dir = "data/s1/s1_sounds/",
              picture_dir = "data/s1/s1_pics/",
              table = df,
              output_dir = "data/s1/",
              output_file = "stimuli_viewer")
```

# Reading data from different linguistic sources

##  Read linguistic files into R

* `textgrid_to_df()` (Praat)
* `eaf_to_df()` (ELAN)
* `exb_to_df()` (EXMARaLDA)
* `srt_to_df()` (subtitle file)
* `audacity_to_df()` (Audacity)
* `flextext_to_df()` (FieldWorks)

##  Read linguistic files into R

```{r}
draw_sound(file_name = "data/test.wav",
           annotation = eaf_to_df("data/test.eaf"))
```


# Get help and cite

## Get help and cite

You can always write an email or open an [issue on GitHub](https://github.com/agricolamz/phonfieldwork/issues), asking some questions.

The most recent citation information is avalible with this command:

```{r}
citation("phonfieldwork")
```

# References {.allowframebreaks}
---
title: "Ethical Research with `phonfieldwork`"
author: "George Moroz, [NRU HSE Linguistic Convergence Laboratory](https://ilcl.hse.ru/en/)"
date: "`r Sys.Date()`"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Ethical Research with `phonfieldwork`}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
```

## Intro

This document is a comment on the ethical policy for the `phonfieldwork` package that does not have any legal effect.

The `create_viewer()` function that is part of the `phonfieldwork` package creates an `.html` file that contains:

* a table with data
* a list of sounds that can be played by users
* a list of pictures that can be viewed by users

The main goal of the `phonfieldwork` package is to create a tool for linguistic research, so I would like to emphasize the possible ethical problems connected to the possibility of putting the obtained `.html` file online. These problems could also be of concern to the ethical comissions of your institution and publishing platforms. Perhapse they have the same or a similar set of rules and concerns as listed below.

## Linguistic research

If you collected data from human subjects I expect that all your participants (or their legal representatives):

* participated voluntarily
* knew the goal of the research
* were informed about possibility to withdraw the consent
* agreed to the obtained data being used for the purpose of scientific research
* agreed to the obtained data being shared with other researchers
* agreed to to the obtained data being made available online (optional, but important if you want to post data online)
* agreed to the proposed compensation

This kind of information is usually collected via informed consent forms, where you also specify the form of data that will be researched and shared: raw data (e.g. audio), aggregated (by speaker, gender, age, etc), anonymized/non-anonymized etc.

So the ethical use of the `phonfieldwork` package implies two things:

* the researcher will not reveal any data that are not listed in the informed consent form 
* the informed consent form does not presuppose the collection of:
    * information related to vulnerable populations that can bring about possible harm
    * personally identifiable or sensitive data that could be traced back to the owner
    
If the subjects in your research do not consent to the publication of their data online but agree to sharing it among researchers, you can use the `encrypt_html_file()` function from the [`encryptedRmd`](https://CRAN.R-project.org/package=encryptedRmd) package in order to make your work password protected.

It is also important to have a contact information in your `.html` in case your subjects will want to withdraw they consent on sharing data.

## Bioacoustic research

It is possible to use phonfieldwork in biacoustic research, so publishing of high quality recordings should be done with caution. Unfortunately, animal call recordings could be used by poachers, hunters, photographers or any other animal lovers in order to lure animals. The possible harm of poachers and hunters is obviouse. Stimulation of animals made by photographers or animal lovers could cause lots of stress and distract from animals' routines (e. g. eating, sleeping etc.). Ethical use of the `phonfieldwork` package implies that the researcher will not reveal any information related to animal populations that may bring about possible harm to any individual animal.

Thanks to Ines Moran for providing information about ethics and bioacoustics.

## Other ethical problems

Since `phonfieldwork` can create an `.html` with any type of data in it, it is important to emphasize other ethical problems not conected to any kind of research. You should not publish data that:

* were produced or distributed without informed consent
* contain, promote, or threaten harassment or violence against other people on the basis of race, ethnicity, national origin, caste, sexual orientation, gender identity, religious affiliation, age, disability, or disease
* contain, or promote suicide or self-harm
* contain, promote or threaten sexual exploitation
* contain, or promote consumption of illegal goods or services
* are synthetic or manipulated in order to cause harm
* violate others’ intellectual property rights, including copyright and trademark
---
title: "Manipulating `phonfieldwork` data with `tidyverse`"
author: "George Moroz, [NRU HSE Linguistic Convergence Laboratory](https://ilcl.hse.ru/en/)"
date: "2021-03-24"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Manipulating `phonfieldwork` data with `tidyverse`}
  %\VignetteEncoding{UTF-8}
  %\VignetteEngine{knitr::rmarkdown}
editor_options: 
  chunk_output_type: console
---



## Introduction

The output tables from the `phonfieldwork`'s functions (e. g. `textgrid_to_df()`, `eaf_to_df()`, `exb_to_df()`, `flextext_to_df()` etc.) is hard to interpret since each row is a separate morpheme/observation or aother unit. In order to merge it to a more common representation we will use `tidyverse` functions (mainly `dplyr` and `tidyr` packages). This text will help you to achieve some results, but it is better to spend some times learning data manipulation with `dplyr` and `tidyr`.

If you do not have `tidyverse` installed run:


```r
install.packages("tidyverse")
```

Let's load the package:


```r
library("tidyverse")
```


## .TextGrid, .eaf, .exb formats
The standard sound annotation formats consisnt of tiers with parallel annotation:

![](unnamed-chunk-31-1.png)

If we convert this file to R we will achieve something like this:


```r
textgrid_to_df("s1/s1_all.TextGrid")
#>    id time_start  time_end      content tier     tier_name          source
#> 1   1  0.0000000 0.4821542          tip    1        labels s1_all.TextGrid
#> 4   1  0.0000000 0.4821542 1_s1_tip.wav    2 backup labels s1_all.TextGrid
#> 7   1  0.0000000 0.1072426                 3               s1_all.TextGrid
#> 8   2  0.1072426 0.1887230            ı    3               s1_all.TextGrid
#> 9   3  0.1887230 0.4821542                 3               s1_all.TextGrid
#> 2   2  0.4821542 0.9120635          tap    1        labels s1_all.TextGrid
#> 5   2  0.4821542 0.9120635 2_s1_tap.wav    2 backup labels s1_all.TextGrid
#> 10  4  0.4821542 0.5770552                 3               s1_all.TextGrid
#> 11  5  0.5770552 0.6793392            æ    3               s1_all.TextGrid
#> 12  6  0.6793392 0.9120635                 3               s1_all.TextGrid
#> 3   3  0.9120635 1.3942177          top    1        labels s1_all.TextGrid
#> 6   3  0.9120635 1.3942177 3_s1_top.wav    2 backup labels s1_all.TextGrid
#> 13  7  0.9120635 1.0364661                 3               s1_all.TextGrid
#> 14  8  1.0364661 1.1066780            ɒ    3               s1_all.TextGrid
#> 15  9  1.1066780 1.3942177                 3               s1_all.TextGrid
```

As we see this table has a long format structure: each observation has its own row. We can select the first two rows with the `filter()` function, remove all unnecessary columns with the `select()` function and spread everything in a table with the `pivot_wider()` function:


```r
textgrid_to_df("s1/s1_all.TextGrid") %>% 
  filter(tier %in% 1:2) %>% 
  select(-time_start, -time_end, -tier_name) %>% 
  pivot_wider(names_from = tier, values_from = content)
#> # A tibble: 3 x 4
#>      id source          `1`   `2`         
#>   <dbl> <chr>           <chr> <chr>       
#> 1     1 s1_all.TextGrid tip   1_s1_tip.wav
#> 2     2 s1_all.TextGrid tap   2_s1_tap.wav
#> 3     3 s1_all.TextGrid top   3_s1_top.wav
```

## .flextext format

Imagine that we obtained the first result from `flextext_to_df()`:


```r
df <- flextext_to_df("files/zilo_test.flextext")
head(df)
#>   p_id s_id w_id    txt     cf hn     gls                   msa                free_trans
#> 1    1    1    1     б-     б-  1      an Inflects any category Жил-был (у Гъули?) петух.
#> 2    1    1    1    ик1    ик1  1    быть                    гл Жил-был (у Гъули?) петух.
#> 3    1    1    1     -о     -о  1     pst               гл:Past Жил-был (у Гъули?) петух.
#> 4    1    1    1     -й     -й  5 cvb(pf)    гл:Converb/Perfect Жил-был (у Гъули?) петух.
#> 5    1    1    1 =гъоди =гъоди  1    =rep                  част Жил-был (у Гъули?) петух.
#> 6    1    1    2     б-     б-  1      an Inflects any category Жил-был (у Гъули?) петух.
#>                            text_title                                morph
#> 1 2017.04 Fairytale about the rooster d7f713db-e8cf-11d3-9764-00c04f186933
#> 2 2017.04 Fairytale about the rooster d7f713e8-e8cf-11d3-9764-00c04f186933
#> 3 2017.04 Fairytale about the rooster d7f713dd-e8cf-11d3-9764-00c04f186933
#> 4 2017.04 Fairytale about the rooster d7f713dd-e8cf-11d3-9764-00c04f186933
#> 5 2017.04 Fairytale about the rooster d7f713e1-e8cf-11d3-9764-00c04f186933
#> 6 2017.04 Fairytale about the rooster d7f713db-e8cf-11d3-9764-00c04f186933
#>                                   word                               phrase
#> 1 efafb420-e203-4685-9be2-1b7810f10a70 1cbadc4f-4051-4783-a0d8-bfeee2d2fb13
#> 2 efafb420-e203-4685-9be2-1b7810f10a70 1cbadc4f-4051-4783-a0d8-bfeee2d2fb13
#> 3 efafb420-e203-4685-9be2-1b7810f10a70 1cbadc4f-4051-4783-a0d8-bfeee2d2fb13
#> 4 efafb420-e203-4685-9be2-1b7810f10a70 1cbadc4f-4051-4783-a0d8-bfeee2d2fb13
#> 5 efafb420-e203-4685-9be2-1b7810f10a70 1cbadc4f-4051-4783-a0d8-bfeee2d2fb13
#> 6 c76d26b7-b84a-42a8-ba34-38e712b1db13 1cbadc4f-4051-4783-a0d8-bfeee2d2fb13
#>                              paragraph                                 text
#> 1 0c9ffe63-b4bf-4af3-a1da-f68567e03513 f08dd466-fca6-4597-925c-c46309387ef7
#> 2 0c9ffe63-b4bf-4af3-a1da-f68567e03513 f08dd466-fca6-4597-925c-c46309387ef7
#> 3 0c9ffe63-b4bf-4af3-a1da-f68567e03513 f08dd466-fca6-4597-925c-c46309387ef7
#> 4 0c9ffe63-b4bf-4af3-a1da-f68567e03513 f08dd466-fca6-4597-925c-c46309387ef7
#> 5 0c9ffe63-b4bf-4af3-a1da-f68567e03513 f08dd466-fca6-4597-925c-c46309387ef7
#> 6 0c9ffe63-b4bf-4af3-a1da-f68567e03513 f08dd466-fca6-4597-925c-c46309387ef7
```

As we can see from `df` print there are three indices in the dataset: `p_id` -- paragraph id, `s_id` -- sentence id and `w_id` -- word id.


```r
df %>% 
  filter(free_trans != "") %>% 
  select(p_id, s_id, w_id, txt, gls, free_trans) %>% 
  group_by(p_id, s_id, free_trans, w_id) %>% 
  summarise(txt = str_c(txt, collapse = ""),
         gls = str_c(gls, collapse = "-"))
#> `summarise()` has grouped output by 'p_id', 's_id', 'free_trans'. You can override using the `.groups` argument.
#> # A tibble: 136 x 6
#> # Groups:   p_id, s_id, free_trans [19]
#>     p_id  s_id free_trans                        w_id txt              gls                       
#>    <dbl> <dbl> <chr>                            <dbl> <chr>            <chr>                     
#>  1     1     1 Жил-был (у Гъули?) петух.            1 б-ик1-о-й=гъоди  "an-быть-pst-cvb(pf)-=rep"
#>  2     1     1 Жил-был (у Гъули?) петух.            2 б--о-ч1игу=гъоди "an--pst-neg.cvb-=rep"    
#>  3     1     1 Жил-был (у Гъули?) петух.            3 Гъули-б          "Гъули-an(gen)"           
#>  4     1     1 Жил-был (у Гъули?) петух.            4 х1елеко          "петух"                   
#>  5     1     1 Жил-был (у Гъули?) петух.            5 .                ""                        
#>  6     2     2 Он грелся на улице(?).               6 къват1и-ла=гъоди "улица-in-=rep"           
#>  7     2     2 Он грелся на улице(?).               7 б-ик1-о-j        "an-быть-pst-cvb(pf)"     
#>  8     2     2 Он грелся на улице(?).               8 букьир-ъа        "Букир-sup"               
#>  9     2     2 Он грелся на улице(?).               9 .                ""                        
#> 10     2     3 [Ему в ногу] воткнулась колючка.    10 къинни-й=гъоди   "втыкаться-cvb(pf)-=rep"  
#> # … with 126 more rows
```

The first `filter()` removes some garbage rows that are present in our example flextext. The `select()` function selects only six important columns from 15 presented in the dataset. The `group_by()` and `summarise()` merge all text from `txt` variable and all glosses from `gls` variable together. Pipe operater `%>% ` make it possible to pass the result from the previous funstion as an input to the following one.

So now we can use the same code in order to merge everything into sentences:


```r
df %>% 
  filter(free_trans != "") %>% 
  select(p_id, s_id, w_id, txt, gls, free_trans) %>% 
  group_by(p_id, s_id, free_trans, w_id) %>% 
  summarise(txt = str_c(txt, collapse = ""),
         gls = str_c(gls, collapse = "-")) %>% 
  group_by(p_id, s_id, free_trans) %>% 
  summarise(txt = str_c(txt, collapse = " "),
         gls = str_c(gls, collapse = " "))
#> `summarise()` has grouped output by 'p_id', 's_id', 'free_trans'. You can override using the `.groups` argument.
#> `summarise()` has grouped output by 'p_id', 's_id'. You can override using the `.groups` argument.
#> # A tibble: 19 x 5
#> # Groups:   p_id, s_id [19]
#>     p_id  s_id free_trans                     txt                        gls                       
#>    <dbl> <dbl> <chr>                          <chr>                      <chr>                     
#>  1     1     1 "Жил-был (у Гъули?) петух."    б-ик1-о-й=гъоди б--о-ч1иг… "an-быть-pst-cvb(pf)-=rep…
#>  2     2     2 "Он грелся на улице(?)."       къват1и-ла=гъоди б-ик1-о-… "улица-in-=rep an-быть-ps…
#>  3     2     3 "[Ему в ногу] воткнулась колю… къинни-й=гъоди ццана .     "втыкаться-cvb(pf)-=rep к…
#>  4     3     4 "Когда колючка воткнулась, [о… ццана къинни-рбигьи б-uʔ-… "колючка втыкаться-pst.pt…
#>  5     4     5 "Гъули не обнаружил дома Бихт… бихьтай=ло ишуишу й-ис-он… "Бихтай-=add дома-дома f-…
#>  6     5     6 "Оттуда пошел к Умалаю, "      б-uʔ-oн-ни=гъоди гьербади… "an-идти-pst-pst(aor)-=re…
#>  7     6     7 "Оттуда петух пошел к Патимат… х1елеко гье-лъу-кку б-ел1… "петух dem-dat-el an-go-p…
#>  8     8    10 "Оттуда [петух] пошел к Ханич… гье-лъу-кку б--и-й б-uʔ-o… "dem-dat-el an--pst(aor)-…
#>  9     9    11 "Иди к Хурмат, ..."            хъаничай-ди=ло ен-л1и беж… "Ханичай-erg-=add rfl-gen…
#> 10    10    12 "Когда дошёл до двора Хурмат,… рул1-и-й Х1урмати-л1и рей… "говорить--cvb(pf) Хурмат…
#> 11    11    13 "Три дня не ели, мы с ней не … рул1-и-й лъоб-гу зубу в-у… "говорить--cvb(pf) три-nu…
#> 12    12    14 "Оттуда он ушёл и дошёл до Ай… гье-лъу-кку б--и-йб--и-й … "dem-dat-el an--pst(aor)-…
#> 13    13    15 "Захраил …?"                   й--и-й Загьраъил-ди=ло б-… "f--pst(aor)-cvb(pf) Захр…
#> 14    14    16 "И он пошел в село. Захраил с… б-укъ-и-й гьеге-б=ло гьон… "an-гнать-pst-pf dem-an-=…
#> 15    15    17 "Оттуда снизу вверх к Исрапил… гьегелъу-кку гьикьу=ло лъ… "там-el внизу-=add вверх …
#> 16    16    18 "Шли-шли и пришли к Гаджи."    гье-лъу-кку б--и-йб--и-й … "dem-dat-el an--pst(aor)-…
#> 17    17    19 "Они поссорились (?) и прогна… й-ейхъ-у й-ах-о-й дунял б… "f-ругать-pst f-драться-p…
#> 18    18    20 "Когда закончили ссориться, [… джид-ия сабаб=ло б-ул1-и-… "делать-fut причина-=add …
#> 19    19    21 "На воротах Забита петух обна… х1елеко б--и-й забити-б к… "петух an--pst(aor)-cvb(p…
```

It is also very easy to get some simple statistics from the data:


```r
df %>% 
  filter(gls != "") %>% 
  count(gls) %>% 
  top_n(6)
#> Selecting by n
#>        gls  n
#> 1     =add 46
#> 2     ¬an1 34
#> 3       an 49
#> 4  cvb(pf) 63
#> 5      pst 28
#> 6 pst(aor) 74
```

Here with the `filter()` function we remove all empty glosses, then we calculate and sort them according to frequency, and in the end we take top six glosses with the `top_n()` function. We can even visualis it with the `ggplot2` package:


```r
df %>% 
  filter(gls != "") %>% 
  count(gls) %>% 
  top_n(11) %>% 
  ggplot(aes(n, fct_reorder(gls, n)))+
  geom_col()+
  labs(x = "count", y = "gloss")
#> Selecting by n
```

![](unnamed-chunk-10-1.png)

---
title: "Phonetic fieldwork and experiments with `phonfieldwork` package"
author: "G. Moroz, [NRU HSE Linguistic Convergence Laboratory](https://ilcl.hse.ru/en/)"
bibliography: bibliography.bib
date: "2021-03-24"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Phonetic fieldwork and experiments with `phonfieldwork` package}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---



# Introduction

There are a lot of different typical tasks that have to be solved during phonetic research and experiments. This includes creating a presentation that will contain all stimuli, renaming and concatenating multiple sound files recorded during a session, automatic annotation in 'Praat' TextGrids (this is one of the sound annotation standards provided by 'Praat' software, see Boersma & Weenink 2020 <https://www.fon.hum.uva.nl/praat/>), creating an html table with annotations and spectrograms, and converting multiple formats ('Praat' TextGrid, 'EXMARaLDA' @schmidt09 and 'ELAN' @wittenburg06). All of these tasks can be solved by a mixture of different tools (any programming language has programs for automatic renaming, and Praat contains scripts for concatenating and renaming files, etc.). `phonfieldwork` provides a functionality that will make it easier to solve those tasks independently of any additional tools. You can also compare the functionality with other packages: ['rPraat'](https://CRAN.R-project.org/package=rPraat) @boril16, ['textgRid'](https://CRAN.R-project.org/package=textgRid) @reidy16, ['pympi'](https://dopefishh.github.io/pympi/index.html) @lubbers13 (thx to Lera Dushkina and Anya Klezovich for letting me know about `pympi`).

There are a lot of different books about linguistic fieldwork and experiments (e.g. @gordon03, @bowern15). This tutorial covers only the data organization part. I will focus on cases where the researcher clearly knows what she or he wants to analyze and has already created a list of stimuli that she or he wants to record. For now `phonfieldwork` works only with `.wav(e)` and `.mp3` audiofiles and `.TextGrid`, `.eaf`, `.exb`, `.srt`, Audacity `.txt` and `.flextext` annotation formats, but the main functionality is availible for `.TextGrid` files (I plan to extend its functionality to other types of data). In the following sections I will describe my workflow for phonetic fieldwork and experiments.

# Install the package

Before you start, make sure that you have installed the package, for example with the following command:


```r
install.packages("phonfieldwork")
```

This command will install the last stable version of the `phonfieldwork` package from CRAN. Since CRAN runs multiple package checks before making it available, this is the safest option. Alternatively, you can download the development version from GitHub:


```r
install.packages("remotes")
remotes::install_github("ropensci/phonfieldwork")
```

If you have any trouble installing the package, you will not be able to use its functionality. In that case you can create [an issue on Github](https://github.com/ropensci/phonfieldwork/issues) or send an email. Since this package could completely destroy your data, **please do not use it until you are sure that you have made a backup**.

Use the `library()` command to load the package:

```r
library("phonfieldwork")
```
In order to work with some `rmarkdown` functions you will need to install `pandoc`, see `vignette("pandoc")` for the details.

This tutorial was made using the following version of `phonfieldwork`:

```r
packageVersion("phonfieldwork")
```

```
## [1] '0.0.12'
```

This tutorial can be cited as follows:

```r
citation("phonfieldwork")
```

```
## 
## Moroz G (2020). _Phonetic fieldwork and experiments with phonfieldwork package_. <URL:
## https://CRAN.R-project.org/package=phonfieldwork>.
## 
## A BibTeX entry for LaTeX users is
## 
##   @Manual{,
##     title = {Phonetic fieldwork and experiments with phonfieldwork package},
##     author = {George Moroz},
##     year = {2020},
##     url = {https://CRAN.R-project.org/package=phonfieldwork},
##   }
```

If you have any trouble using the package, do not hesitate to create [an issue on Github](https://github.com/ropensci/phonfieldwork/issues/new).


# Philosophy of the `phonfieldwork` package

Most phonetic research consists of the following steps:

1. Formulate a research question. Think of what kind of data is necessary to answer this question, what is the appropriate amount of data, what kind of annotation you will do, what kind of statistical models and visualizations you will use, etc.
2. Create a list of stimuli.
3. Elicite list of stimuli from speakers who signed an *Informed Consent* statement, agreeing to participate in the experiment to be recorded on audio and/or video. Keep an eye on recording settings: sampling rate, resolution (bit), and number of channels should be the same across all recordings.
4. Annotate the collected data.
5. Extract the collected data.
6. Create visualizations and evaluate your statistical models.
7. Report your results.
8. Publish your data.

The `phonfieldwork` package is created for helping with items 3, partially with 4, and 5 and 8.

To make the automatic annotation of data easier, I usually record each stimulus as a separate file.  While recording, I carefully listen to my consultants to make sure that they are producing the kind of speech I want: three isolated pronunciations of the same stimulus, separated by a pause and contained in a carrier phrase. In case a speaker does not produce three clear repetitions, I ask them to repeat the task, so that as a result of my fieldwork session I will have:

* a collection of small soundfiles (video) with the same sampling rate, resolution (bit), and number of channels
* a list of succesful and unsuccesful attempts to produce a stimulus according to my requirements (usually I keep this list in a regular notebook)

There are some phoneticians who prefer to record everything, for language documentation purposes. I think that should be a separate task: you can’t have your cake and eat it too. But if you insist on recording everything, it is possible to run two recorders at the same time: one could run during the whole session, while the other is used to produce small audio files. You can also use special software to record your stimuli automatically on a computer (e.g. [SpeechRecorder](https://www.bas.uni-muenchen.de/Bas/software/speechrecorder/) or [PsychoPy](https://www.psychopy.org/)).

You can show a native speaker your stimuli one by one or not show them the stimule but ask them to pronounce a certain stimulus or its translation. I use presentations to collect all stimuli in a particular order without the risk of omissions.

Since each stimulus is recorded as a separate audiofile, it is possible to merge them into one file automatically and make an annotation in a Praat TextGrid (the same result can be achieved with the `Concatenate recoverably` command in Praat). After this step, the user needs to do some annotation of her/his own. When the annotation part is finished, it is possible to extract the annotated parts to a table, where each annotated object is a row characterised by some features (stimulus, repetition, speaker, etc...). You can play the soundfile and view its oscilogram and spectrogram. Here is [an example of such a file](https://ropensci.github.io/phonfieldwork/additional/stimuli_viewer.html) and [instruction for doing it](https://ropensci.github.io/phonfieldwork/articles/phonfieldwork.html#create-a-viewer-1).

# The `phonfieldwork` package in use
## Make a list of your stimuli

There are several ways to enter information about a list of stimuli into R:

* using the `c()` function you can create a **vector** of all words and store it in a variable `my_stimuli` (you can choose any other name):


```r
my_stimuli <- c("tip", "tap", "top")
```

* it is also possible to store your list as a column in a `.csv` file and read it into R using the `read.csv()` function:


```r
my_stimuli_df <- read.csv("my_stimuli_df.csv")
my_stimuli_df
```

```
##   stimuli vowel
## 1     tip     ı
## 2     tap     æ
## 3     top     ɒ
```

* it is also possible to store your list as a column in an `.xls` or `xlsx` file and read it into R using the `read_xls` or `read_xlsx` functions from the `readxl` package. If the package `readxl` is not installed on your computer, install it using `install.packages("readxl")`


```r
library("readxl")
# run install.packages("readxl") in case you don't have it installed
my_stimuli_df <- read_xlsx("my_stimuli_df.xlsx")
my_stimuli_df
```

```
## # A tibble: 3 x 2
##   stimuli vowel
##   <chr>   <chr>
## 1 tip     ı    
## 2 tap     æ    
## 3 top     ɒ
```

## Create a presentation based on a list of stimuli

When the list of stimuli is loaded into R, you can create a presentation for elicitation. It is important to define an output directory, so in the following example I use the `getwd()` function, which returns the path to the current working directory. You can set any directory as your current one using the `setwd()` function. It is also possible to provide a path to your intended output directory with `output_dir` (e. g. "/home/user_name/..."). This command (unlike `setwd()`) does not change your working directory.


```r
create_presentation(stimuli = my_stimuli_df$stimuli,
                    output_file = "first_example",
                    output_dir = getwd())
```

As a result, a file "first_example.html" was created in the output folder. You can change the name of this file by changing  the `output_file` argument. The `.html` file now looks as follows:

<iframe src="https://ropensci.github.io/phonfieldwork/additional/first_example.html" width = 800 height = 600>
  <p>Your browser does not support iframes :(</p>
</iframe>

It is also possible to change the output format, using the `output_format` argument. By dafault it is "html", but you can also use "pptx" (this is a relatively new feature of `rmarkdown`, so update the package in case you get errors). There is also an additional argument `translations`, where you can provide translations for stimuli in order that they appeared near the stimuli on the slide.

It is also possible to use images (or gif, e. g. for a sign language research) as a stimuli. In order to do that you need to provide an absolute or relative path to the file instead of the stimulus and mark in the `external` argument, which of the stimuli is external:


```r
my_image <- system.file("extdata", "r-logo.png", package = "phonfieldwork")
my_image
```

```
## [1] "/home/agricolamz/R/x86_64-pc-linux-gnu-library/4.0/phonfieldwork/extdata/r-logo.png"
```

```r
create_presentation(stimuli = c("rzeka", "drzewo", my_image),
                    external = 3,
                    output_file = "second_example",
                    output_dir = getwd())
```

<iframe src="https://ropensci.github.io/phonfieldwork/additional/second_example.html" width = 800 height = 600>
  <p>Your browser does not support iframes :(</p>
</iframe>

## Rename collected data
After collecting data and removing soundfiles with unsuccesful elicitations, one could end up with the following structure:


```
## ├── s1
## │   ├── 01.wav
## │   ├── 02.wav
## │   └── 03.wav
## ├── s2
## │   ├── 01.wav
## │   ├── 02.wav
## │   └── 03.wav
```

For each speaker `s1` and `s2` there is a folder that containes three audiofiles. Now let's rename the files.


```r
rename_soundfiles(stimuli = my_stimuli_df$stimuli,
                  prefix = "s1_",
                  path = "s1/")
```

```
## You can find change correspondences in the following file:
## /home/agricolamz/work/packages/phonfieldwork/vignettes/s1/backup/logging.csv
```

As a result, you obtain the following structure:


```
## ├── s1
## │   ├── 1_s1_tip.wav
## │   ├── 2_s1_tap.wav
## │   ├── 3_s1_top.wav
## │   └── backup
## │       ├── 01.wav
## │       ├── 02.wav
## │       ├── 03.wav
## │       └── logging.csv
## ├── s2
## │   ├── 01.wav
## │   ├── 02.wav
## │   └── 03.wav
```

The `rename_soundfiles()` function created a backup folder with all of the unrenamed files, and renamed all files using the prefix provided in the `prefix` argument. There is an additional argument `backup` that can be set to `FALSE` (it is `TRUE` by default), in case you are sure that the renaming function will work properly with your files and stimuli, and you do not need a backup of the unrenamed files. There is also an additional argument `logging` (`TRUE` by default) that creates a `logging.csv` file in the `backup` folder (or in the original folder if the `backup` argument has value `FALSE`) with the correspondences between old and new names of the files. Here is the contence of the `logging.csv`:


```
##     from           to
## 1 01.wav 1_s1_tip.wav
## 2 02.wav 2_s1_tap.wav
## 3 03.wav 3_s1_top.wav
```

To each name was added an additional prefix with number that make it easear to keep the original sorting of the stimuli. If you do not want this autonumbering turn the `autonumbering` to `FALSE`:


```r
rename_soundfiles(stimuli = my_stimuli_df$stimuli,
                  prefix = "s2_",
                  suffix = paste0("_", 1:3),
                  path = "s2/",
                  backup = FALSE,
                  logging = FALSE,
                  autonumbering = FALSE)
```


```
## ├── s1
## │   ├── 1_s1_tip.wav
## │   ├── 2_s1_tap.wav
## │   ├── 3_s1_top.wav
## │   └── backup
## │       ├── 01.wav
## │       ├── 02.wav
## │       ├── 03.wav
## │       └── logging.csv
## ├── s2
## │   ├── s2_tap_2.wav
## │   ├── s2_tip_1.wav
## │   └── s2_top_3.wav
```

The last command renamed the soundfiles in the `s2` folder, adding the prefix `s2` as in the previous example, and the suffix `1`-`3`. On most operating systems it is impossible to create two files with the same name, so sometimes it can be useful to add some kind of index at the end of the files.

There is also a possible scenario, that not all stimuli are retrieved from informant. So in order to deal with that case there is an additional argument `missing`, where user can put id numbers of stimuli that are not present in audiofiles:


```r
rename_soundfiles(stimuli = my_stimuli_df$stimuli,
                  path = "s3/",
                  missing = c(1, 3))
```

Sometimes it is useful to get information about sound duration:


```r
get_sound_duration("s1/2_s1_tap.wav")
```

```
##           file  duration
## 1 2_s1_tap.wav 0.4821542
```

It is also possible to analyze the whole folder using the `read_from_folder()` function. The first argument is the path to the folder. The second argument is the type of information or file type (possible values: "audacity", "duration", "eaf", "exb", "flextext", "formant", "intensity", "picth", "srt", "textgrid"):


```r
read_from_folder(path = "s2/", "duration")
```

```
##           file  duration
## 1 s2_tap_2.wav 0.5343991
## 2 s2_tip_1.wav 0.5866440
## 3 s2_top_3.wav 0.6650113
```

For now `phonfieldwork` works only with `.wav(e)` and `.mp3` sound files.

## Merge all data together

After all the files are renamed, you can merge them into one. Remmber that sampling rate, resolution (bit), and number of channels should be the same across all recordings. It is possible to resample files with the `resample()` function from  `biacoustics`.


```r
concatenate_soundfiles(path = "s1/",
                       result_file_name = "s1_all")
```

This comand creates a new soundfile `s1_all.wav` and an asociated Praat TextGrid `s1_all.TextGrid`:


```
## ├── s1
## │   ├── 1_s1_tip.wav
## │   ├── 2_s1_tap.wav
## │   ├── 3_s1_top.wav
## │   ├── backup
## │   │   ├── 01.wav
## │   │   ├── 02.wav
## │   │   ├── 03.wav
## │   │   └── logging.csv
## │   ├── s1_all.TextGrid
## │   └── s1_all.wav
## ├── s2
## │   ├── s2_tap_2.wav
## │   ├── s2_tip_1.wav
## │   └── s2_top_3.wav
```

The resulting file can be parsed with Praat:

![](unnamed-chunk-23-1.png)

Sometimes recorded sounds do not have any silence at the beginning or the end, so after the merging the result utterances will too close to each other. It is possible to fix using the argument `separate_duration` of the `concatenate_soundfiles()` function: just put the desired duration of the separator in seconds.

It is not kind of task that could occur within phonfieldwork philosophy, but it also possible to merge multiple `.TextGrid`s with the same tier structure using `concatente_textgrids()` function.

## Annotate your data

It is possible to annotate words using an existing annotation:


```r
my_stimuli_df$stimuli
```

```
## [1] "tip" "tap" "top"
```

```r
annotate_textgrid(annotation =  my_stimuli_df$stimuli,
                  textgrid = "s1/s1_all.TextGrid")
```

![](unnamed-chunk-25-1.png)

As you can see in the example, the `annotate_textgrid()` function creates a backup of the tier and adds a new tier on top of the previous one. It is possible to prevent the function from doing so by setting the `backup` argument to `FALSE`.

Imagine that we are interested in annotation of vowels. The most common solution will be open Praat and create new annotations. But it is also possible to create them in advance using subannotations. The idea that you choose some baseline tier that later will be automatically cutted into smaller pieces on the other tier.


```r
create_subannotation(textgrid = "s1/s1_all.TextGrid",
                     tier = 1, # this is a baseline tier
                     n_of_annotations = 3) # how many empty annotations per unit?
```

![](unnamed-chunk-27-1.png)

It is worth mentioning that if you want to have different number of subannotation per unit, you can pass a vector of required numbers to `n_of_annotations` argument.

After the creation of subannotations, we can annotate created tier:


```r
annotate_textgrid(annotation = c("", "ı", "", "", "æ", "", "", "ɒ", ""),
                  textgrid = "s1/s1_all.TextGrid",
                  tier = 3,
                  backup = FALSE)
```

![](unnamed-chunk-29-1.png)

You can see that we created a third tier with annotation. The only thing left is to move annotation boundaries in Praat (this can not be automated):



![](unnamed-chunk-31-1.png)

You can see from the last figure that no backup tier was created (`backup = FALSE`), that the third tier was annotated (`tier = 3`).

In case you want to create an empty TextGrid it is possible to use a `create_empty_textgrid()` function that takes a duration as an argument:


```r
create_empty_textgrid(get_sound_duration("s2/s2_tip_1.wav")$duration,
                      tier_name = c("a", "b"),
                      path = "s2",
                      result_file_name = "s2_tip_1")
```
![](unnamed-chunk-33-1.png)

```
## ├── s1
## │   ├── 1_s1_tip.wav
## │   ├── 2_s1_tap.wav
## │   ├── 3_s1_top.wav
## │   ├── backup
## │   │   ├── 01.wav
## │   │   ├── 02.wav
## │   │   ├── 03.wav
## │   │   └── logging.csv
## │   ├── s1_all.TextGrid
## │   └── s1_all.wav
## ├── s2
## │   ├── s2_tap_2.wav
## │   ├── s2_tip_1.TextGrid
## │   ├── s2_tip_1.wav
## │   └── s2_top_3.wav
```

It is also possible to remove some tier from textgrid. For instance, we can remove one tier from the previously created file:


```r
remove_textgrid_tier(textgrid = "s2/s2_tip_1.TextGrid", tier = 2)
```
![](unnamed-chunk-36-1.png)

## Extracting your data

First, it is important to create a folder where all of the extracted files will be stored:

```r
dir.create("s1/s1_sounds")
```

It is possible to extract all annotated files based on an annotation tier:


```r
extract_intervals(file_name = "s1/s1_all.wav",
                  textgrid = "s1/s1_all.TextGrid",
                  tier = 3,
                  path = "s1/s1_sounds/",
                  prefix = "s1_")
```


```
## ├── s1
## │   ├── 1_s1_tip.wav
## │   ├── 2_s1_tap.wav
## │   ├── 3_s1_top.wav
## │   ├── backup
## │   │   ├── 01.wav
## │   │   ├── 02.wav
## │   │   ├── 03.wav
## │   │   └── logging.csv
## │   ├── s1_all.TextGrid
## │   ├── s1_all.wav
## │   └── s1_sounds
## │       ├── 1_s1_ı.wav
## │       ├── 2_s1_æ.wav
## │       └── 3_s1_ɒ.wav
## ├── s2
## │   ├── s2_tap_2.wav
## │   ├── s2_tip_1.TextGrid
## │   ├── s2_tip_1.wav
## │   └── s2_top_3.wav
```

## Visualizing your data
It is possible to view an oscilogram and spetrogram of any soundfile:


```r
draw_sound(file_name = "s1/s1_sounds/1_s1_ı.wav")
```

![](unnamed-chunk-40-1.png)

There are additional parameters:

* `title` -- the title for the plot
* `from` -- time in seconds at which to start extraction
* `to` -- time in seconds at which to stop extraction
* `zoom` -- time in seconds for zooming spectrogram
* `text_size` -- size of the text on the plot
* `annotation` -- the optional file with the TextGrid's file path or dataframe with annotations (see the section 5.)
* `freq_scale` -- the measure of the frequency: can be "Hz" or "kHz".
* `frequency_range` -- the frequency range to be displayed for the spectrogram
* `dynamic_range` -- values greater than this many dB below the maximum will be displayed in the same color
* `window_length` -- the desired length in milliseconds for the analysis window
* `window` -- window type (can be "rectangular", "hann", "hamming", "cosine", "bartlett", "gaussian", and "kaiser")
* `preemphasisf` -- Preemphasis of 6 dB per octave is added to frequencies above the specified frequency. For no preemphasis (important for bioacoustics), set to a 0.
* `spectrum_info` -- logical value, if `FALSE` won't print information about spectorgram on the right side of the plot.
* `output_file` -- the name of the output file
* `output_width` -- the width of the device
* `output_height` -- the height of the device
* `output_units` -- the units in which height and width are given. This can be "px" (pixels, which is the default value), "in" (inches), "cm" or "mm".

It is really important in case you have a long file not to draw the whole file, since it won't fit into the RAM of your computer. So you can use `from` and `to` arguments in order to plot the fragment of the sound and annotation:


```r
draw_sound("s1/s1_all.wav",
           "s1/s1_all.TextGrid",
           from = 0.4,
           to = 0.95)
```

```
## Warning in df$tier == unique(df$tier): longer object length is not a multiple of shorter object
## length
```

![](unnamed-chunk-41-1.png)

It is also possible using the `zoom` argument to show the part of the spectrogram keeping the broader oscilogram context:


```r
draw_sound("s1/s1_all.wav",
           "s1/s1_all.TextGrid",
           zoom = c(0.4, 0.95))
```

```
## Warning in df$tier == unique(df$tier): longer object length is not a multiple of shorter object
## length
```

![](unnamed-chunk-42-1.png)

If the `output_file` argument is provided, R will save the plot in your directory instead of displaying it.


```r
draw_sound(file_name = "s1/s1_sounds/1_s1_ı.wav",
           output_file = "s1/s1_tip",
           title = "s1 tip")
```


```
## ├── s1
## │   ├── 1_s1_tip.wav
## │   ├── 2_s1_tap.wav
## │   ├── 3_s1_top.wav
## │   ├── backup
## │   │   ├── 01.wav
## │   │   ├── 02.wav
## │   │   ├── 03.wav
## │   │   └── logging.csv
## │   ├── s1_all.TextGrid
## │   ├── s1_all.wav
## │   ├── s1_sounds
## │   │   ├── 1_s1_ı.wav
## │   │   ├── 2_s1_æ.wav
## │   │   └── 3_s1_ɒ.wav
## │   └── s1_tip.png
## ├── s2
## │   ├── s2_tap_2.wav
## │   ├── s2_tip_1.TextGrid
## │   ├── s2_tip_1.wav
## │   └── s2_top_3.wav
```

It is also possible to create visualizations of all sound files in a folder. For this purpose you need to specify a source folder with the argument `sounds_from_folder` and a target folder for the images (`pic_folder_name`). The new image folder is automatically created in the upper level folder, so that sound and image folders are on the same level in the tree structure of your directory.


```r
draw_sound(sounds_from_folder = "s1/s1_sounds/",
           pic_folder_name = "s1_pics")
```


```
## ├── s1
## │   ├── 1_s1_tip.wav
## │   ├── 2_s1_tap.wav
## │   ├── 3_s1_top.wav
## │   ├── backup
## │   │   ├── 01.wav
## │   │   ├── 02.wav
## │   │   ├── 03.wav
## │   │   └── logging.csv
## │   ├── s1_all.TextGrid
## │   ├── s1_all.wav
## │   ├── s1_pics
## │   │   ├── 1_s1_ı.png
## │   │   ├── 2_s1_æ.png
## │   │   └── 3_s1_ɒ.png
## │   ├── s1_sounds
## │   │   ├── 1_s1_ı.wav
## │   │   ├── 2_s1_æ.wav
## │   │   └── 3_s1_ɒ.wav
## │   └── s1_tip.png
## ├── s2
## │   ├── s2_tap_2.wav
## │   ├── s2_tip_1.TextGrid
## │   ├── s2_tip_1.wav
## │   └── s2_top_3.wav
```

It is also possible to use the argument `textgrid_from_folder` in order to specify the folder where .TextGrids for annotation are (could be the same folder as the sound one). By default the `draw_sound()` function with the `sounds_from_folder` argument adds a title with the file name to each pictures' title, but it is possible to turn it off using the argument `title_as_filename = FALSE`.

If you are familiar with the Raven program for bioacoustics, you probably miss an ability to annotate not only time, but also a frequency range. In order to do it you need to create a dataframe with the columns `time_start`, `time_end`, `freq_low` and `freq_high`:


```r
raven_an <- data.frame(time_start = 450,
                       time_end  = 520,
                       freq_low = 3,
                       freq_high = 4.5)

draw_sound(system.file("extdata", "test.wav", package = "phonfieldwork"),
           raven_annotation = raven_an)
```

![](unnamed-chunk-47-1.png)

It is also possible to use multiple values, colors (adding `colors` column) and annotation (adding `content` column):


```r
raven_an <- data.frame(time_start = c(250, 450),
                       time_end  = c(400, 520),
                       freq_low = c(1, 3),
                       freq_high = c(2, 4.5),
                       colors = c("red", "blue"),
                       content = c("a", "b"))

draw_sound(system.file("extdata", "test.wav", package = "phonfieldwork"),
           raven_annotation = raven_an)
```

![](unnamed-chunk-48-1.png)

# Read linguistic files into R

The `phonfieldwork` package provides also several methods for reading different file types into R. This makes it possible to analyze them and convert into `.csv` files (e. g. using the `write.csv()` function). The main advantage of using those functions is that all of them return `data.frame`s with columns (`time_start`, `time_end`, `content` and `source`). This make it easer to use the result in the `draw_sound()` function that make it possible to visualise all kind of sound annotation systems.

* file `.TextGrid` from Praat (just change the `system.file()` function to path to the file); see also [`rPraat`](https://fu.ff.cuni.cz/praat/#rpraat-package-for-r) and [`textgRid`](https://github.com/patrickreidy/textgRid) packages

```r
textgrid_to_df(system.file("extdata", "test.TextGrid", package = "phonfieldwork"))
```

```
##    id time_start   time_end content tier       tier_name        source
## 1   1 0.00000000 0.01246583            1       intervals test.TextGrid
## 6   1 0.00000000 0.01246583            2 empty_intervals test.TextGrid
## 2   2 0.01246583 0.24781914       t    1       intervals test.TextGrid
## 7   2 0.01246583 0.24781914            2 empty_intervals test.TextGrid
## 11  1 0.01246583 0.01246583       t    3          points test.TextGrid
## 3   3 0.24781914 0.39552363       e    1       intervals test.TextGrid
## 8   3 0.24781914 0.39552363            2 empty_intervals test.TextGrid
## 12  2 0.24781914 0.24781914       e    3          points test.TextGrid
## 4   4 0.39552363 0.51157715       s    1       intervals test.TextGrid
## 9   4 0.39552363 0.51157715            2 empty_intervals test.TextGrid
## 13  3 0.39552363 0.39552363       s    3          points test.TextGrid
## 5   5 0.51157715 0.65267574       t    1       intervals test.TextGrid
## 10  5 0.51157715 0.65267574            2 empty_intervals test.TextGrid
## 14  4 0.51157715 0.51157715       t    3          points test.TextGrid
```

* file `.eaf` from ELAN  (just change the `system.file()` function to path to the file); see also the [FRelan](https://github.com/langdoc/FRelan) package by Niko Partanen


```r
eaf_to_df(system.file("extdata", "test.eaf", package = "phonfieldwork"))
```

```
##    tier id content       tier_name tier_type time_start time_end   source
## 9     1  1               intervals     praat      0.000    0.012 test.eaf
## 10    2  1         empty_intervals     praat      0.000    0.012 test.eaf
## 11    1  2       t       intervals     praat      0.012    0.248 test.eaf
## 12    2  2       C empty_intervals     praat      0.012    0.248 test.eaf
## 1     1  3       e       intervals     praat      0.248    0.396 test.eaf
## 2     2  3       V empty_intervals     praat      0.248    0.396 test.eaf
## 3     1  4       s       intervals     praat      0.396    0.512 test.eaf
## 4     2  4       C empty_intervals     praat      0.396    0.512 test.eaf
## 5     1  5       t       intervals     praat      0.512    0.652 test.eaf
## 6     2  5       C empty_intervals     praat      0.512    0.652 test.eaf
## 7     1  6               intervals     praat      0.652  300.000 test.eaf
## 8     2  6         empty_intervals     praat      0.652  300.000 test.eaf
```

* file `.exb` from EXMARaLDA  (just change the `system.file()` function to path to the file)

```r
exb_to_df(system.file("extdata", "test.exb", package = "phonfieldwork"))
```

```
##   tier id content tier_name tier_type tier_category tier_speaker time_start  time_end   source
## 3    1  1       t     X [v]         t             v         SPK0 0.06908955 0.2498984 test.exb
## 1    1  2       e     X [v]         t             v         SPK0 0.24989836 0.3807275 test.exb
## 5    1  3       s     X [v]         t             v         SPK0 0.38072750 0.4042473 test.exb
## 7    1  4       t     X [v]         t             v         SPK0 0.40424735 0.6526757 test.exb
## 4    2  1       C     X [v]         a             v         SPK0 0.06908955 0.2498984 test.exb
## 2    2  2       V     X [v]         a             v         SPK0 0.24989836 0.3807275 test.exb
## 6    2  3       C     X [v]         a             v         SPK0 0.38072750 0.4042473 test.exb
## 8    2  4       C     X [v]         a             v         SPK0 0.40424735 0.6526757 test.exb
```

* subtitles file `.srt` (just change the `system.file()` function to path to the file)


```r
srt_to_df(system.file("extdata", "test.srt", package = "phonfieldwork"))
```

```
##   id content time_start time_end   source
## 0  1       t      0.013    0.248 test.srt
## 1  2       e      0.248    0.396 test.srt
## 2  3       s      0.396    0.512 test.srt
## 3  4       t      0.512    0.653 test.srt
```


* file `.txt` from Audacity


```r
audacity_to_df(system.file("extdata", "test_audacity.txt", package = "phonfieldwork"))
```

```
##   time_start  time_end content            source
## 1  0.2319977 0.3953891    sssw test_audacity.txt
```


* file `.flextext` from FLEx (that is actually is not connected with the main functionality of `phonfieldwork`, but I'd like to have it):


```r
head(flextext_to_df("files/zilo_test.flextext"))
```

```
##   p_id s_id w_id    txt     cf hn     gls                   msa                free_trans
## 1    1    1    1     б-     б-  1      an Inflects any category Жил-был (у Гъули?) петух.
## 2    1    1    1    ик1    ик1  1    быть                    гл Жил-был (у Гъули?) петух.
## 3    1    1    1     -о     -о  1     pst               гл:Past Жил-был (у Гъули?) петух.
## 4    1    1    1     -й     -й  5 cvb(pf)    гл:Converb/Perfect Жил-был (у Гъули?) петух.
## 5    1    1    1 =гъоди =гъоди  1    =rep                  част Жил-был (у Гъули?) петух.
## 6    1    1    2     б-     б-  1      an Inflects any category Жил-был (у Гъули?) петух.
##                            text_title                                morph
## 1 2017.04 Fairytale about the rooster d7f713db-e8cf-11d3-9764-00c04f186933
## 2 2017.04 Fairytale about the rooster d7f713e8-e8cf-11d3-9764-00c04f186933
## 3 2017.04 Fairytale about the rooster d7f713dd-e8cf-11d3-9764-00c04f186933
## 4 2017.04 Fairytale about the rooster d7f713dd-e8cf-11d3-9764-00c04f186933
## 5 2017.04 Fairytale about the rooster d7f713e1-e8cf-11d3-9764-00c04f186933
## 6 2017.04 Fairytale about the rooster d7f713db-e8cf-11d3-9764-00c04f186933
##                                   word                               phrase
## 1 efafb420-e203-4685-9be2-1b7810f10a70 1cbadc4f-4051-4783-a0d8-bfeee2d2fb13
## 2 efafb420-e203-4685-9be2-1b7810f10a70 1cbadc4f-4051-4783-a0d8-bfeee2d2fb13
## 3 efafb420-e203-4685-9be2-1b7810f10a70 1cbadc4f-4051-4783-a0d8-bfeee2d2fb13
## 4 efafb420-e203-4685-9be2-1b7810f10a70 1cbadc4f-4051-4783-a0d8-bfeee2d2fb13
## 5 efafb420-e203-4685-9be2-1b7810f10a70 1cbadc4f-4051-4783-a0d8-bfeee2d2fb13
## 6 c76d26b7-b84a-42a8-ba34-38e712b1db13 1cbadc4f-4051-4783-a0d8-bfeee2d2fb13
##                              paragraph                                 text
## 1 0c9ffe63-b4bf-4af3-a1da-f68567e03513 f08dd466-fca6-4597-925c-c46309387ef7
## 2 0c9ffe63-b4bf-4af3-a1da-f68567e03513 f08dd466-fca6-4597-925c-c46309387ef7
## 3 0c9ffe63-b4bf-4af3-a1da-f68567e03513 f08dd466-fca6-4597-925c-c46309387ef7
## 4 0c9ffe63-b4bf-4af3-a1da-f68567e03513 f08dd466-fca6-4597-925c-c46309387ef7
## 5 0c9ffe63-b4bf-4af3-a1da-f68567e03513 f08dd466-fca6-4597-925c-c46309387ef7
## 6 0c9ffe63-b4bf-4af3-a1da-f68567e03513 f08dd466-fca6-4597-925c-c46309387ef7
```

There is also an additional function for working with the `.flextext` format that convert it to a glossed document in a `docx` or `.html` format (see examples: [`.docx`](https://github.com/ropensci/phonfieldwork/raw/master/docs/additional/glossed_document.docx), [`.html`](https://ropensci.github.io/phonfieldwork/additional/glossed_document.html)):


```r
create_glossed_document(flextext = "files/zilo_test.flextext",
                        output_dir = ".") # you need to specify the path to the output folder
```

```
## Output created: /home/agricolamz/work/packages/phonfieldwork/vignettes/glossed_document.html
```



It is also possible to convert to LaTeX examples format using the `example_pkg` argument (possible values are: `gb4e`, `langsci`, `expex`, `philex`). There is also an additional [text](https://ropensci.github.io/phonfieldwork/articles/data_manipulation_with_tidyverse.html) about manipulation with `flextext_to_df()` output.

All those functions (`tier_to_df()`, `textgrid_to_df()`, `eaf_to_df()`, `exb_to_df()`, `audacity_to_df()`, `srt_to_df()`) except `flextext_to_df()` can be used in order to visualise sound annotation:


```r
draw_sound(file_name = system.file("extdata", "test.wav", package = "phonfieldwork"),
           annotation = eaf_to_df(system.file("extdata", "test.eaf", package = "phonfieldwork")))
```

![](unnamed-chunk-57-1.png)

Remmember that it is also possible to read multiple files with the `read_from_folder()` funtion.

# Create a viewer

Sound viewer (here is an [example 1](https://ropensci.github.io/phonfieldwork/additional/stimuli_viewer.html) and [example 2](https://ropensci.github.io/phonfieldwork/additional/stimuli_viewer2.html)) is a useful tool that combine together your annotations and make it searchable. It is also produce a ready to go `.html` file that could be uploaded on the server (e. g. to Github Pages) and be availible for anyone in the world.

In order to create a sound viewer you need three things:

* folder with sounds
* folder with pictures
* dataframe with some details (e. g. annotation, utterance number etc.)

We will start with the previous folder structure:


```
## ├── s1
## │   ├── 1_s1_tip.wav
## │   ├── 2_s1_tap.wav
## │   ├── 3_s1_top.wav
## │   ├── backup
## │   │   ├── 01.wav
## │   │   ├── 02.wav
## │   │   ├── 03.wav
## │   │   └── logging.csv
## │   ├── s1_all.TextGrid
## │   ├── s1_all.wav
## │   ├── s1_pics
## │   │   ├── 1_s1_ı.png
## │   │   ├── 2_s1_æ.png
## │   │   └── 3_s1_ɒ.png
## │   ├── s1_sounds
## │   │   ├── 1_s1_ı.wav
## │   │   ├── 2_s1_æ.wav
## │   │   └── 3_s1_ɒ.wav
## │   └── s1_tip.png
## ├── s2
## │   ├── s2_tap_2.wav
## │   ├── s2_tip_1.TextGrid
## │   ├── s2_tip_1.wav
## │   └── s2_top_3.wav
```

We have all folders:


```r
list.files("s1/s1_sounds/") # sounds
```

```
## [1] "1_s1_ı.wav" "2_s1_æ.wav" "3_s1_ɒ.wav"
```

```r
list.files("s1/s1_pics/") # pictures
```

```
## [1] "1_s1_ı.png" "2_s1_æ.png" "3_s1_ɒ.png"
```

So what is left is the table. It is possible to create manually (or upload it form .csv or .xlsx files, see section 4.1):


```r
df <- data.frame(word  = c("tap", "tip", "top"),
                 sounds = c("æ", "ı", "ɒ"))
df
```

```
##   word sounds
## 1  tap      æ
## 2  tip      ı
## 3  top      ɒ
```

This table could be used in order to create an annotation viewer:


```r
create_viewer(audio_dir = "s1/s1_sounds/",
              picture_dir = "s1/s1_pics/",
              table = df,
              output_dir = "s1/",
              output_file = "stimuli_viewer")
```

```
## Since the result .html file possibly containes some vulnerable data, researcher(s) bear the whole responsibility for the publishing of the result. Run vignette("ethical_research_with_phonfieldwork") for more details.
```

```
## Output created: s1/stimuli_viewer.html
```

As a result, a `stimuli_viewer.html` was created in the `s1` folder.


```
## ├── s1
## │   ├── 1_s1_tip.wav
## │   ├── 2_s1_tap.wav
## │   ├── 3_s1_top.wav
## │   ├── backup
## │   │   ├── 01.wav
## │   │   ├── 02.wav
## │   │   ├── 03.wav
## │   │   └── logging.csv
## │   ├── s1_all.TextGrid
## │   ├── s1_all.wav
## │   ├── s1_pics
## │   │   ├── 1_s1_ı.png
## │   │   ├── 2_s1_æ.png
## │   │   └── 3_s1_ɒ.png
## │   ├── s1_sounds
## │   │   ├── 1_s1_ı.wav
## │   │   ├── 2_s1_æ.wav
## │   │   └── 3_s1_ɒ.wav
## │   ├── s1_tip.png
## │   └── stimuli_viewer.html
## ├── s2
## │   ├── s2_tap_2.wav
## │   ├── s2_tip_1.TextGrid
## │   ├── s2_tip_1.wav
## │   └── s2_top_3.wav
```

You can find the created example [here](https://ropensci.github.io/phonfieldwork/additional/stimuli_viewer.html).

Unfortunately, the way of table creation for the annotation viewer presented in this section is not a good solution for the huge amount of sounds. It is possible to derive such a table from annotation TextGrid, that we have created earlier. Here is a TextGrid:


```r
textgrid_to_df("s1/s1_all.TextGrid")
```

```
##    id time_start  time_end      content tier     tier_name          source
## 1   1  0.0000000 0.4821542          tip    1        labels s1_all.TextGrid
## 4   1  0.0000000 0.4821542 1_s1_tip.wav    2 backup labels s1_all.TextGrid
## 7   1  0.0000000 0.1072426                 3               s1_all.TextGrid
## 8   2  0.1072426 0.1887230            ı    3               s1_all.TextGrid
## 9   3  0.1887230 0.4821542                 3               s1_all.TextGrid
## 2   2  0.4821542 0.9120635          tap    1        labels s1_all.TextGrid
## 5   2  0.4821542 0.9120635 2_s1_tap.wav    2 backup labels s1_all.TextGrid
## 10  4  0.4821542 0.5770552                 3               s1_all.TextGrid
## 11  5  0.5770552 0.6793392            æ    3               s1_all.TextGrid
## 12  6  0.6793392 0.9120635                 3               s1_all.TextGrid
## 3   3  0.9120635 1.3942177          top    1        labels s1_all.TextGrid
## 6   3  0.9120635 1.3942177 3_s1_top.wav    2 backup labels s1_all.TextGrid
## 13  7  0.9120635 1.0364661                 3               s1_all.TextGrid
## 14  8  1.0364661 1.1066780            ɒ    3               s1_all.TextGrid
## 15  9  1.1066780 1.3942177                 3               s1_all.TextGrid
```

So in order to create desired table we can use `tier_to_df()` function:


```r
t1 <- tier_to_df("s1/s1_all.TextGrid", tier = 1)
t1
```

```
##   id time_start  time_end content tier tier_name          source
## 1  1  0.0000000 0.4821542     tip    1    labels s1_all.TextGrid
## 2  2  0.4821542 0.9120635     tap    1    labels s1_all.TextGrid
## 3  3  0.9120635 1.3942177     top    1    labels s1_all.TextGrid
```

```r
t3 <- tier_to_df("s1/s1_all.TextGrid", tier = 3)
t3
```

```
##    id time_start  time_end content tier tier_name          source
## 7   1  0.0000000 0.1072426            3           s1_all.TextGrid
## 8   2  0.1072426 0.1887230       ı    3           s1_all.TextGrid
## 9   3  0.1887230 0.4821542            3           s1_all.TextGrid
## 10  4  0.4821542 0.5770552            3           s1_all.TextGrid
## 11  5  0.5770552 0.6793392       æ    3           s1_all.TextGrid
## 12  6  0.6793392 0.9120635            3           s1_all.TextGrid
## 13  7  0.9120635 1.0364661            3           s1_all.TextGrid
## 14  8  1.0364661 1.1066780       ɒ    3           s1_all.TextGrid
## 15  9  1.1066780 1.3942177            3           s1_all.TextGrid
```

As we see the first tier is ready, but the third tier contains empty annotations. Let's remove them:


```r
t3 <- t3[t3$content != "",]
t3
```

```
##    id time_start  time_end content tier tier_name          source
## 8   2  0.1072426 0.1887230       ı    3           s1_all.TextGrid
## 11  5  0.5770552 0.6793392       æ    3           s1_all.TextGrid
## 14  8  1.0364661 1.1066780       ɒ    3           s1_all.TextGrid
```

So from this point it is possible to create the table that we wanted:


```r
new_df <- data.frame(words = t1$content,
                     sounds = t3$content)
new_df
```

```
##   words sounds
## 1   tip      ı
## 2   tap      æ
## 3   top      ɒ
```

So now we are ready to run our code for creating an annotation viewer:


```r
create_viewer(audio_dir = "s1/s1_sounds/",
              picture_dir = "s1/s1_pics/",
              table = new_df,
              output_dir = "s1/",
              output_file = "stimuli_viewer")
```

```
## Since the result .html file possibly containes some vulnerable data, researcher(s) bear the whole responsibility for the publishing of the result. Run vignette("ethical_research_with_phonfieldwork") for more details.
```

```
## Output created: s1/stimuli_viewer.html
```

By default sorting in the result annotation viewer will be according file names in the system, so if you want to have another default sorting you can specify column names that the result table should be sorted by using the `sorting_columns` argument.

If you are familiar with my package [`lingtypology`](https://ropensci.github.io/lingtypology/) @moroz17 for interactive linguistic map generation and API for typological databases, there is a good news for you: it is possible to connect those two pacakages creating an interactive map that share the same hear and view buttons. In order to do it you need

* to add a `glottocode` column with language glottocodes from [Glottolog](https://glottolog.org/) @hammarstrom2020 to your dataframe with annotation details;
* install `lingtypology` with a command `install.packages("lingtypology")` if you don't have it installed;
* add `map = TRUE` argument to `create_viewer()` function.

I will add some glottocodes for Russian, Polish and Czech to the dataframe that we have already worked with (for those data it doesn't make any sense, I just giving an example of usage):


```r
new_df$glottocode <- c("russ1263", "poli1260", "czec1258")
create_viewer(audio_dir = "s1/s1_sounds/",
              picture_dir = "s1/s1_pics/",
              table = new_df,
              output_dir = "s1/",
              output_file = "stimuli_viewer2",
              map = TRUE)
```

```
## Since the result .html file possibly containes some vulnerable data, researcher(s) bear the whole responsibility for the publishing of the result. Run vignette("ethical_research_with_phonfieldwork") for more details.
```

```
## Output created: s1/stimuli_viewer2.html
```

[Here](https://ropensci.github.io/phonfieldwork/additional/stimuli_viewer2.html) is the result file.

It is also possible to provide your own coordinates with `latitude` and `longitude` columns. In that case `glottocode` column is optional.

# References
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/exb_to_df.R
\name{exb_to_df}
\alias{exb_to_df}
\title{EXMARaLDA's .exb file to dataframe}
\usage{
exb_to_df(file_name)
}
\arguments{
\item{file_name}{string with a filename or path to the .exb file}
}
\value{
a dataframe with columns:  \code{tier}, \code{id}, \code{content},
\code{tier_name}, \code{tier_type}, \code{tier_category},
\code{tier_speaker}, \code{time_start}, \code{time_end}, \code{source}.
}
\description{
Convert .exb file from EXMARaLDA to a dataframe.
}
\examples{
exb_to_df(system.file("extdata", "test.exb", package = "phonfieldwork"))
}
\author{
George Moroz <agricolamz@gmail.com>
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/draw_spectrogram.R
\name{draw_spectrogram}
\alias{draw_spectrogram}
\title{Draw spectrograms}
\usage{
draw_spectrogram(
  sound,
  fs = 22050,
  text_size = 1,
  window_length = 5,
  dynamic_range = 50,
  window = "kaiser",
  windowparameter = -1,
  freq_scale = "kHz",
  spectrum_info = TRUE,
  timestep = -1000,
  padding = 10,
  preemphasisf = 50,
  frequency_range = c(0, 5),
  nlevels = dynamic_range,
  x_axis = TRUE,
  title = NULL,
  raven_annotation = NULL,
  formant_df = NULL
)
}
\arguments{
\item{sound}{Either a numeric vector representing a sequence of samples taken
from a sound wave or a sound object created with the loadsound() or
makesound() functions.}

\item{fs}{The sampling frequency in Hz. If a sound object is passed this
does not need to be specified.}

\item{text_size}{numeric, text size (default = 1).}

\item{window_length}{The desired analysis window length in milliseconds.}

\item{dynamic_range}{Values greater than this many dB below the maximum will
be displayed in the same color.}

\item{window}{A string indicating the type of window desired. Supported types
are: rectangular, hann, hamming, cosine, bartlett, gaussian, and kaiser.}

\item{windowparameter}{The parameter necessary to generate the window, if
appropriate. At the moment, the only windows that require parameters are the
Kaiser and Gaussian windows. By default, these are set to 2 for kaiser and
0.4 for gaussian windows.}

\item{freq_scale}{a string indicating the type of frequency scale. Supported
types are: "Hz" and "kHz".}

\item{spectrum_info}{logical. If \code{TRUE} then add information about
window method and params.}

\item{timestep}{If a negative value is given, -N, then N equally-spaced time
steps are calculated. If a positive number is given, this is the spacing
between adjacent analyses, in milliseconds.}

\item{padding}{The amount of zero padding for each window, measured in units
of window length. For example, if the window is 50 points, and padding = 10,
500 zeros will be appended to each window.}

\item{preemphasisf}{Preemphasis of 6 dB per octave is added to frequencies
above the specified frequency. For no preemphasis, set to a frequency higher
than the sampling frequency.}

\item{frequency_range}{vector with the range of frequencies to be displayed
for the spectrogram up to a maximum of \code{fs}/2. This is set to 0-5 kHz by
default.}

\item{nlevels}{The number of divisions to be used for the z-axis of the
spectrogram. By default it is set equal to the dynamic range, meaning that a
single color represents 1 dB on the z-axis.}

\item{x_axis}{If \code{TRUE} then draw x axis.}

\item{title}{Character with the title.}

\item{raven_annotation}{Raven (Center for Conservation Bioacoustics) style
annotations (boxes over spectrogram). The dataframe that contains
\code{time_start}, \code{time_end}, \code{freq_low} and \code{freq_high}
columns. Optional columns are \code{colors} and \code{content}.}

\item{formant_df}{dataframe with formants from \code{formant_to_df()} function}
}
\description{
This function was slightly changed from \code{phonTools::spectrogram()}.
Argument description is copied from \code{phonTools::spectrogram()}.
}
\examples{
\dontrun{
draw_spectrogram(system.file("extdata", "test.wav",
  package = "phonfieldwork"
))
}

}
\author{
Santiago Barreda <sbarreda@ucdavis.edu>
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/read_from_folder.R
\name{read_from_folder}
\alias{read_from_folder}
\title{Read multiple files from the folder}
\usage{
read_from_folder(path, type)
}
\arguments{
\item{path}{to a folder with multiple sound files.}

\item{type}{should be one of the following: "duration", "audacity", "eaf", "exb", "flextext", "formant", "intensity", "picth", "srt", "textgrid"}
}
\description{
This function reads multiple files from the folder. The first argument is the path, the second argument is the type of files to read.
}
\examples{

read_from_folder(system.file("extdata", package = "phonfieldwork"), "eaf")

}
\author{
George Moroz <agricolamz@gmail.com>
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rename_soundfiles.R
\name{rename_soundfiles}
\alias{rename_soundfiles}
\title{Rename soundfiles}
\usage{
rename_soundfiles(
  stimuli,
  translations = NULL,
  prefix = NULL,
  suffix = NULL,
  order = NULL,
  missing = NULL,
  path,
  autonumbering = TRUE,
  backup = TRUE,
  logging = TRUE
)
}
\arguments{
\item{stimuli}{character vector of stimuli}

\item{translations}{character vector of translations (optonal). This values are added after stimuli to the new files' names so the result will be \code{...stimulus_translation...}.}

\item{prefix}{character vector of length one containing prefix for file names}

\item{suffix}{character vector of length one containing suffix for file names}

\item{order}{numeric vector that define the order of stimuli. By default the
order of the stimuli is taken.}

\item{missing}{numeric vector that define missing stimuli in case when some stimuli are not recorded.}

\item{path}{path to the directory with soundfiles.}

\item{autonumbering}{logical. If TRUE, function creates an automatic numbering of files.}

\item{backup}{logical. If TRUE, function creates backup folder with all
files. By default is TRUE.}

\item{logging}{logical. If TRUE creates a .csv file with the correspondences of old names and new names. This could be useful for restoring in case something goes wrong.}
}
\value{
no output
}
\description{
Rename soundfiles using the template from user.
}
\author{
George Moroz <agricolamz@gmail.com>
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/eaf_to_df.R
\name{eaf_to_df}
\alias{eaf_to_df}
\title{ELAN's .eaf file to dataframe}
\usage{
eaf_to_df(file_name)
}
\arguments{
\item{file_name}{string with a filename or path to the .eaf file}
}
\value{
a dataframe with columns:  \code{tier}, \code{id}, \code{content},
\code{tier_name}, \code{tier_type}, \code{time_start}, \code{time_end},
\code{source}).
}
\description{
Convert .eaf file from ELAN to a dataframe.
}
\examples{
eaf_to_df(system.file("extdata", "test.eaf", package = "phonfieldwork"))
}
\author{
George Moroz <agricolamz@gmail.com>
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/create_viewer.R
\name{create_viewer}
\alias{create_viewer}
\title{Create an annotation viewer}
\usage{
create_viewer(
  audio_dir,
  picture_dir,
  table,
  captions = NULL,
  sorting_columns = NULL,
  about = "Created with the `phonfieldworks` package (Moroz 2020).",
  map = FALSE,
  output_dir,
  output_file = "stimuli_viewer",
  render = TRUE
)
}
\arguments{
\item{audio_dir}{path to the directory with sounds}

\item{picture_dir}{path to the directory with pictures}

\item{table}{data frame with data ordered according to files in the audio folder}

\item{captions}{vector of strings that will be used for captions for a picture.}

\item{sorting_columns}{vector of strings for sorting the result column}

\item{about}{it is either .Rmd file or string with the text for about information: author, project, place of gahtered information and other metadata, version of the viewer and so on}

\item{map}{the logical argument, if \code{TRUE} and there is a \code{glottocode} column in \code{table}}

\item{output_dir}{the output directory for the rendered file}

\item{output_file}{the name of the result .html file (by default stimuli_viewer)}

\item{render}{the logical argument, if \code{TRUE} renders the created R Markdown viewer to the \code{output_dir} folder, otherwise returns the path to the temporary file with a .csv file.}
}
\value{
If \code{render} is \code{FALSE}, the function returns a path to the temporary file with .csv file. If \code{render} is \code{TRUE}, there is no output in a function.
}
\description{
Creates an html file with table and sound preview and player
}
\author{
George Moroz <agricolamz@gmail.com>
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_textgrid_names.R
\name{get_textgrid_names}
\alias{get_textgrid_names}
\title{Extract TextGrid names}
\usage{
get_textgrid_names(textgrid)
}
\arguments{
\item{textgrid}{path to the TextGrid}
}
\value{
return a vector of tier names from given TextGrid
}
\description{
Extract TextGrid names.
}
\examples{
get_textgrid_names(system.file("extdata", "test.TextGrid",
  package = "phonfieldwork"
))
}
\author{
George Moroz <agricolamz@gmail.com>
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tier_to_df.R
\name{tier_to_df}
\alias{tier_to_df}
\title{TextGrid's tier to dataframe}
\usage{
tier_to_df(file_name, tier = 1)
}
\arguments{
\item{file_name}{string with a filename or path to the TextGrid}

\item{tier}{value that could be either ordinal number of the tier either
name of the tier. By default is '1'.}
}
\value{
a dataframe with columns:  \code{id}, \code{time_start},
\code{time_end}, \code{content}, , \code{tier_name}
}
\description{
Convert selected tier from a Praat TextGrid to a dataframe.
}
\examples{
tier_to_df(system.file("extdata", "test.TextGrid",
  package = "phonfieldwork"
))
tier_to_df(
  system.file("extdata", "test.TextGrid",
    package = "phonfieldwork"
  ),
  "intervals"
)
}
\author{
George Moroz <agricolamz@gmail.com>
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/create_empty_textgrid.R
\name{create_empty_textgrid}
\alias{create_empty_textgrid}
\title{Create an empty TextGrid}
\usage{
create_empty_textgrid(
  duration,
  tier_name = NULL,
  point_tier = NULL,
  path,
  result_file_name = "new_textgrid"
)
}
\arguments{
\item{duration}{integer. Duration of the textgrid. If you do not know the duration of your audio file use the \code{get_sound_duration()} function.}

\item{tier_name}{a vector that contain tier names.}

\item{point_tier}{a vector that defines which tiers should be made point tiers. This argument excepts numeric values (e. g. \code{c(2, 4)} means second and forth tiers) or character (e. g. \code{c("a", "b")} means tiers with names "a" and "b")}

\item{path}{path to the directory with soundfiles.}

\item{result_file_name}{name of the result and annotation files.}
}
\value{
The function returns no output, just creates a Praat TextGrid in the same folder as a reference sound file.
}
\description{
Creates an empty Praat TextGrid in the same folder as a reference sound file. It is possible to manage with predefined number of tiers, their names and their types.
}
\examples{
tmp <- tempfile(fileext = ".TextGrid")
create_empty_textgrid(1, path = dirname(tmp), result_file_name = basename(tmp))

}
\author{
George Moroz <agricolamz@gmail.com>
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/formant_to_df.R
\name{formant_to_df}
\alias{formant_to_df}
\title{Praat Formant object to dataframe}
\usage{
formant_to_df(file_name)
}
\arguments{
\item{file_name}{string with a filename or path to the Formant file}
}
\value{
a dataframe with columns:  \code{time_start}, \code{time_end},
\code{frequency}, \code{bandwidth} and \code{formant}
}
\description{
Convert a Praat Formant object to a dataframe.
}
\examples{
formant_to_df(system.file("extdata", "e.Formant", package = "phonfieldwork"))
}
\author{
George Moroz <agricolamz@gmail.com>
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/draw_sound.R
\name{draw_sound}
\alias{draw_sound}
\title{Draw Oscilogram, Spectrogram and annotation}
\usage{
draw_sound(
  file_name,
  annotation = NULL,
  from = NULL,
  to = NULL,
  zoom = NULL,
  text_size = 1,
  output_file = NULL,
  title = NULL,
  freq_scale = "kHz",
  frequency_range = c(0, 5),
  dynamic_range = 50,
  window_length = 5,
  window = "kaiser",
  windowparameter = -1,
  preemphasisf = 50,
  spectrum_info = TRUE,
  raven_annotation = NULL,
  formant_df = NULL,
  pitch = NULL,
  pitch_range = c(75, 350),
  intensity = NULL,
  output_width = 750,
  output_height = 500,
  output_units = "px",
  sounds_from_folder = NULL,
  textgrids_from_folder = NULL,
  pic_folder_name = "pics",
  title_as_filename = TRUE,
  prefix = NULL,
  suffix = NULL,
  autonumber = FALSE
)
}
\arguments{
\item{file_name}{a sound file}

\item{annotation}{a source for annotation files (path to TextGrid file or dataframe created from other linguistic types, e. g. via \code{textgrid_to_df()}, \code{eaf_to_df()} or other functions)}

\item{from}{Time in seconds at which to start extraction.}

\item{to}{Time in seconds at which to stop extraction.}

\item{zoom}{numeric vector of zoom window time (in seconds). It will draw
the whole oscilogram and part of the spectrogram.}

\item{text_size}{numeric, text size (default = 1).}

\item{output_file}{the name of the output file}

\item{title}{the title for the plot}

\item{freq_scale}{a string indicating the type of frequency scale.
Supported types are: "Hz" and "kHz".}

\item{frequency_range}{vector with the range of frequencies to be displayed
for the spectrogram up to a maximum of fs/2. By default this is set to 0-5
kHz.}

\item{dynamic_range}{values greater than this many dB below the maximum will
be displayed in the same color}

\item{window_length}{the desired analysis window length in milliseconds.}

\item{window}{A string indicating the type of window desired. Supported types
are: "rectangular", "hann", "hamming", "cosine", "bartlett", "gaussian", and
"kaiser".}

\item{windowparameter}{The parameter necessary to generate the window, if
appropriate. At the moment, the only windows that require parameters are the
Kaiser and Gaussian windows. By default, these are set to 2 for kaiser and
0.4 for gaussian windows.}

\item{preemphasisf}{Preemphasis of 6 dB per octave is added to frequencies
above the specified frequency. For no preemphasis, set to a frequency higher
than the sampling frequency.}

\item{spectrum_info}{logical. If \code{TRUE} then add information about
window method and params.}

\item{raven_annotation}{Raven (Center for Conservation Bioacoustics) style
annotations (boxes over spectrogram). The dataframe that contains
\code{time_start}, \code{time_end}, \code{freq_low} and \code{freq_high}
columns. Optional columns are \code{colors} and \code{content}.}

\item{formant_df}{dataframe with formants from \code{formant_to_df()} function}

\item{pitch}{path to the Praat `.Pitch` file or result of
\code{pitch_to_df()} function. This variable provide data for visualisation
of a pitch contour exported from Praat.}

\item{pitch_range}{vector with the range of frequencies to be displayed.
By default this is set to 75-350 Hz.}

\item{intensity}{path to the Praat `.Intensity` file or result of
\code{intensity_to_df()} function. This variable provide data for
visualisation of an intensity contour exported from Praat.}

\item{output_width}{the width of the device}

\item{output_height}{the height of the device}

\item{output_units}{the units in which height and width are given.
Can be "px" (pixels, the default), "in" (inches), "cm" or "mm".}

\item{sounds_from_folder}{path to a folder with multiple sound files.
If this argument is not \code{NULL}, then the function goes through all
files and creates picture for all of them.}

\item{textgrids_from_folder}{path to a folder with multiple .TextGrid files.
If this argument is not \code{NULL}, then the function goes through all files
and create picture for all of them.}

\item{pic_folder_name}{name for a folder, where all pictures will be stored
in case \code{sounds_from_folder} argument is not \code{NULL}}

\item{title_as_filename}{logical. If true adds filename title to each picture}

\item{prefix}{prefix for all file names for created pictures in case
\code{sounds_from_folder} argument is not \code{NULL}}

\item{suffix}{suffix for all file names for created pictures in case
\code{sounds_from_folder} argument is not \code{NULL}}

\item{autonumber}{if TRUE automatically add number of extracted sound to the
file_name. Prevents from creating a duplicated files and wrong sorting.}
}
\value{
Oscilogram and spectrogram plot (and possibly TextGrid annotation).
}
\description{
Create oscilogram and spectrogram plot.
}
\examples{
\dontrun{
draw_sound(system.file("extdata", "test.wav", package = "phonfieldwork"))

draw_sound(
  system.file("extdata", "test.wav", package = "phonfieldwork"),
  system.file("extdata", "test.TextGrid",
    package = "phonfieldwork"
  )
)

draw_sound(system.file("extdata", "test.wav", package = "phonfieldwork"),
  system.file("extdata", "test.TextGrid", package = "phonfieldwork"),
  pitch = system.file("extdata", "test.Pitch",
    package = "phonfieldwork"
  ),
  pitch_range = c(50, 200)
)
draw_sound(system.file("extdata", "test.wav", package = "phonfieldwork"),
  system.file("extdata", "test.TextGrid", package = "phonfieldwork"),
  pitch = system.file("extdata", "test.Pitch",
    package = "phonfieldwork"
  ),
  pitch_range = c(50, 200),
  intensity = intensity_to_df(system.file("extdata", "test.Intensity",
    package = "phonfieldwork"
  ))
)
draw_sound(system.file("extdata", "test.wav", package = "phonfieldwork"),
  formant_df = formant_to_df(system.file("extdata", "e.Formant",
    package = "phonfieldwork"
  ))
)
}
}
\author{
George Moroz <agricolamz@gmail.com>
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/intensity_to_df.R
\name{intensity_to_df}
\alias{intensity_to_df}
\title{Praat Intensity tier to dataframe}
\usage{
intensity_to_df(file_name)
}
\arguments{
\item{file_name}{string with a filename or path to the Intensity tier}
}
\value{
a dataframe with columns:  \code{time_start}, \code{time_end},
\code{Intensity}
}
\description{
Convert a Praat Intensity tier to a dataframe.
}
\examples{
intensity_to_df(system.file("extdata", "test.Intensity", package = "phonfieldwork"))
}
\author{
George Moroz <agricolamz@gmail.com>
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/create_presentation.R
\name{create_presentation}
\alias{create_presentation}
\title{Creates a presentation}
\usage{
create_presentation(
  stimuli,
  translations = "",
  external = NULL,
  font_size = 50,
  output_dir,
  output_format = "html",
  output_file = "stimuli_presentation",
  render = TRUE
)
}
\arguments{
\item{stimuli}{the vector of stimuli (obligatory). Can be a path to an image.}

\item{translations}{the vector of translations (optional)}

\item{external}{the vector with the indices of external images}

\item{font_size}{font size in px (50, by default)}

\item{output_dir}{the output directory for the rendered file}

\item{output_format}{the string that difine the R Markdown output format:
"html" (by default) or "pptx"}

\item{output_file}{the name of the result presentation file
(by default stimuli_presentation)}

\item{render}{the logical argument, if \code{TRUE} render the created R
Markdown presentation to the \code{output_dir} folder, otherwise returns the
path to the temporary file with a Rmd file.}
}
\value{
If \code{render} is \code{FALSE}, the function returns a path to the
temporary file. If \code{render} is \code{TRUE}, there is no output in a
function.
}
\description{
Creates an html or powerpoint presentation in a working directory from list
of words and translations.
\href{https://ropensci.github.io/phonfieldwork/additional/first_example.html}{Here}
is an example of such presentation.
}
\examples{
create_presentation(
  stimuli = c("rzeka", "drzewo"),
  translations = c("river", "tree"),
  render = FALSE
)

# with image
create_presentation(
  stimuli = c(
    "rzeka", "drzewo",
    system.file("extdata", "r-logo.png",
      package = "phonfieldwork"
    )
  ),
  translations = c("river", "tree", ""),
  external = 3,
  render = FALSE
)
}
\author{
George Moroz <agricolamz@gmail.com>
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/concatenate_textgrids.R
\name{concatenate_textgrids}
\alias{concatenate_textgrids}
\title{Concatenate sounds}
\usage{
concatenate_textgrids(path, result_file_name = "concatenated")
}
\arguments{
\item{path}{path to the directory with soundfiles.}

\item{result_file_name}{name of the result and annotation files.}
}
\value{
no output
}
\description{
Creates a merged sound file from old sound files in a folder. If the annotation argument is not equal to \code{NULL}, it creates an annotation file (Praat .TextGrid, ELAN .eaf or EXMARaLDA .exb) with original sound names annotation.
}
\examples{
# create two files in a temprary folder "test_folder"
t1 <- system.file("extdata", "test.TextGrid", package = "phonfieldwork")
t2 <- system.file("extdata", "post.TextGrid", package = "phonfieldwork")
tdir <- tempdir()
file.copy(c(t1, t2), tdir)

# here are two .wav files in a folder
list.files(tdir)
# [1] "post.TextGrid" "test.TextGrid" ...

# Concatenate all TextGrids from the folder into concatenated.TextGrid
concatenate_textgrids(path = tdir, result_file_name = "concatenated")

list.files(tdir)
# [1] "concatenated.TextGrid" "post.TextGrid" "test.TextGrid" ...
}
\author{
George Moroz <agricolamz@gmail.com>
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.r
\name{create_image_look_up}
\alias{create_image_look_up}
\title{Create image look_up objects for html viewer}
\usage{
create_image_look_up(img_src, img_caption = NULL, text = "&#x1f441;")
}
\arguments{
\item{img_src}{string or vector of strings with a image(s) path(s).}

\item{img_caption}{string or vector of strings that will be displayed when
image is clicked.}

\item{text}{string o vector of strings that will be displayed as view link.
By default it is eye emoji (&#x1f441;).}
}
\value{
a string or vector of strings
}
\description{
Create image look_up objects for html viewer
}
\author{
George Moroz <agricolamz@gmail.com>
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/annotate_textgrid.R
\name{annotate_textgrid}
\alias{annotate_textgrid}
\title{Annotate textgrid}
\usage{
annotate_textgrid(
  annotation,
  textgrid,
  tier = 1,
  each = 1,
  backup = TRUE,
  write = TRUE
)
}
\arguments{
\item{annotation}{vector of stimuli}

\item{textgrid}{character with a filename or path to the TextGrid}

\item{tier}{value that could be either ordinal number of the tier either name
of the tier}

\item{each}{non-negative integer. Each element of x is repeated each times}

\item{backup}{logical. If TRUE (by default) it creates a backup tier.}

\item{write}{logical. If TRUE (by dafault) it overwrites an existing tier.}
}
\value{
a string that contain TextGrid. If argument write is \code{TRUE},
then no output.
}
\description{
Annotates textgrids. It is possible to define step in the argument "each",
so each second element of the tier will be annotated.
}
\examples{
annotate_textgrid(
  annotation = c("", "t", "e", "s", "t"),
  textgrid = system.file("extdata",
    "test.TextGrid",
    package = "phonfieldwork"
  ),
  tier = 2, write = FALSE
)
}
\author{
George Moroz <agricolamz@gmail.com>
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/flextext_to_df.R
\name{flextext_to_df}
\alias{flextext_to_df}
\title{FLEX's .flextext file to dataframe}
\usage{
flextext_to_df(file_name)
}
\arguments{
\item{file_name}{string with a filename or path to the .flextext file}
}
\value{
a dataframe with columns: \code{p_id}, \code{s_id}, \code{w_id},
\code{txt}, \code{cf}, \code{hn}, \code{gls},
\code{msa}, \code{morph}, \code{word}, \code{phrase}, \code{paragraph},
\code{free_trans}, \code{text}, \code{text_title}
}
\description{
Convert .flextext file from FLEX to a dataframe.
}
\author{
George Moroz <agricolamz@gmail.com>
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.r
\name{create_sound_play}
\alias{create_sound_play}
\title{Create audio play objects for html viewer}
\usage{
create_sound_play(snd_src, text = "&#x1f442;")
}
\arguments{
\item{snd_src}{string or vector of strings with a image(s) path(s).}

\item{text}{string o vector of strings that will be displayed as view link.
By default it is ear emoji (&#x1f442;).}
}
\value{
a string or vector of strings
}
\description{
Create audio play objects for html viewer
}
\author{
George Moroz <agricolamz@gmail.com>
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pitch_to_df.R
\name{pitch_to_df}
\alias{pitch_to_df}
\title{Praat Pitch tier to dataframe}
\usage{
pitch_to_df(file_name, candidates = "")
}
\arguments{
\item{file_name}{string with a filename or path to the Pitch tier}

\item{candidates}{Praat Pitch tier contains multiple candidates for each
time slice, use the value \code{"all"} if you want to get them all}
}
\value{
a dataframe with columns:  \code{time_start}, \code{time_end},
\code{frequency} and, if \code{candidates} = \code{"all"},
\code{candidate_id} and \code{strength}
}
\description{
Convert a Praat Pitch tier to a dataframe.
}
\examples{
pitch_to_df(system.file("extdata", "test.Pitch", package = "phonfieldwork"))
}
\author{
George Moroz <agricolamz@gmail.com>
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/concatenate_soundfiles.R
\name{concatenate_soundfiles}
\alias{concatenate_soundfiles}
\title{Concatenate sounds}
\usage{
concatenate_soundfiles(
  path,
  result_file_name = "concatenated",
  annotation = "textgrid",
  separate_duration = 0
)
}
\arguments{
\item{path}{path to the directory with soundfiles.}

\item{result_file_name}{name of the result and annotation files.}

\item{annotation}{character. There are several variants: "textgrid" for Praat TextGrid, "eaf" for ELAN's .eaf file, or "exb" for EXMARaLDA's .exb file. It is also possible to use \code{NULL} in order to prevent the creation of the annotation file.}

\item{separate_duration}{double. It is possible to add some silence between concatenated sounds. This variable denotes duration of this soundless separator in seconds.}
}
\value{
no output
}
\description{
Creates a merged sound file from old sound files in a folder. If the annotation argument is not equal to \code{NULL}, it creates an annotation file (Praat .TextGrid, ELAN .eaf or EXMARaLDA .exb) with original sound names annotation.
}
\examples{
# create two files in a temprary folder "test_folder"
s1 <- system.file("extdata", "test.wav", package = "phonfieldwork")
s2 <- system.file("extdata", "post.wav", package = "phonfieldwork")
tdir <- tempdir()
file.copy(c(s1, s2), tdir)

# here are two .wav files in a folder
list.files(tdir)
# [1] "post.wav" "test.wav" ...

# Concatenate all files from the folder into concatenated.wav and create
# corresponding TextGrid
concatenate_soundfiles(path = tdir, result_file_name = "concatenated")

list.files(tdir)
# [1] "concatenated.TextGrid" "concatenated.wav" "post.wav" "test.wav" ...
}
\author{
George Moroz <agricolamz@gmail.com>
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/set_textgrid_names.R
\name{set_textgrid_names}
\alias{set_textgrid_names}
\title{Rewrite TextGrid names}
\usage{
set_textgrid_names(textgrid, tiers, names, write = TRUE)
}
\arguments{
\item{textgrid}{path to the TextGrid}

\item{tiers}{integer vector with the number of tiers that should be named}

\item{names}{vector of strings with new names for TextGrid tiers}

\item{write}{logical. If TRUE (by dafault) it overwrites an existing tier}
}
\value{
a string that contain TextGrid. If argument write is \code{TRUE},
then no output.
}
\description{
Rewrite TextGrid names.
}
\examples{
set_textgrid_names(system.file("extdata", "test.TextGrid",
  package = "phonfieldwork"
),
tiers = 3, names = "new_name", write = FALSE
)
}
\author{
George Moroz <agricolamz@gmail.com>
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_sound_duration.R
\name{get_sound_duration}
\alias{get_sound_duration}
\title{Get file(s) duration}
\usage{
get_sound_duration(file_name)
}
\arguments{
\item{file_name}{a sound file}
}
\description{
Calculate sound(s) duration.
}
\examples{
get_sound_duration(
  system.file("extdata", "test.wav", package = "phonfieldwork")
)
}
\author{
George Moroz <agricolamz@gmail.com>
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/create_glossed_document.R
\name{create_glossed_document}
\alias{create_glossed_document}
\title{Create a glossed document}
\usage{
create_glossed_document(
  flextext = NULL,
  rows = c("gls"),
  output_dir,
  output_file = "glossed_document",
  output_format = "html",
  example_pkg = NULL
)
}
\arguments{
\item{flextext}{path to a .flextext file or a dataframe with the following
columns: \code{p_id}, \code{s_id}, \code{w_id}, \code{txt}, \code{cf},
\code{hn}, \code{gls}, \code{msa}, \code{morph}, \code{word}, \code{phrase},
\code{paragraph}, \code{free_trans}, \code{text}, \code{text_title}}

\item{rows}{vector of row names from the flextext that should appear in the
final document. Possible values are: "cf", "hn", "gls", "msa". "gls" is
default.}

\item{output_dir}{the output directory for the rendered file}

\item{output_file}{the name of the result \code{.html} file (by default
\code{glossed_document}).}

\item{output_format}{The option can be "html" or "docx"}

\item{example_pkg}{vector with name of the LaTeX package for glossing
(possible values: \code{"gb4e"}, \code{"langsci"}, \code{"expex"},
\code{"philex"})}
}
\value{
If \code{render} is \code{FALSE}, the function returns a path to
the temporary file with .csv file. If \code{render} is \code{TRUE}, there is
no output in a function.
}
\description{
Creates a file with glossed example (export from .flextext or other formats)
}
\author{
George Moroz <agricolamz@gmail.com>
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/srt_to_df.R
\name{srt_to_df}
\alias{srt_to_df}
\title{Subtitles .srt file to dataframe}
\usage{
srt_to_df(file_name)
}
\arguments{
\item{file_name}{string with a filename or path to the .srt file}
}
\value{
a dataframe with columns:  \code{id}, \code{content},
\code{time_start}, \code{time_end}, \code{source}.
}
\description{
Convert subtitles .srt file to a dataframe.
}
\examples{
srt_to_df(system.file("extdata", "test.srt", package = "phonfieldwork"))
}
\author{
George Moroz <agricolamz@gmail.com>
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/audacity_to_df.R
\name{audacity_to_df}
\alias{audacity_to_df}
\title{Audacity's labels to dataframe}
\usage{
audacity_to_df(file_name)
}
\arguments{
\item{file_name}{file_name string with a filename or path to the .txt file
produced by Audacity}
}
\value{
a dataframe with columns:  \code{content}, \code{time_start},
\code{time_end}, \code{source}.
}
\description{
Audacity make it possible to annotate sound files with labels that can be
exported as a .tsv file with .txt extension. This function convert result to
dataframe.
}
\examples{
audacity_to_df(system.file("extdata",
  "test_audacity.txt",
  package = "phonfieldwork"
))
}
\author{
George Moroz <agricolamz@gmail.com>
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/df_to_tier.R
\name{df_to_tier}
\alias{df_to_tier}
\title{Dataframe to TextGrid's tier}
\usage{
df_to_tier(df, textgrid, tier_name = "", overwrite = TRUE)
}
\arguments{
\item{df}{an R dataframe object that contains columns named "content",
"time_start" and "time_end"}

\item{textgrid}{a character with a filename or path to the TextGrid}

\item{tier_name}{a vector that contain a name for a created tier}

\item{overwrite}{a logic argument, if \code{TRUE} overwrites the existing
TextGrid file}
}
\value{
If \code{overwrite} is \code{FALSE}, then the function returns a
vector of strings with a TextGrid. If \code{overwrite} is \code{TRUE}, then
no output.
}
\description{
Convert a dataframe to a Praat TextGrid.
}
\examples{
time_start <- c(0.00000000, 0.01246583, 0.24781914, 0.39552363, 0.51157715)
time_end <- c(0.01246583, 0.24781914, 0.39552363, 0.51157715, 0.65267574)
content <- c("", "T", "E", "S", "T")
df_to_tier(my_df <- data.frame(id = 1:5, time_start, time_end, content),
  system.file("extdata", "test.TextGrid",
    package = "phonfieldwork"
  ),
  overwrite = FALSE
)
}
\author{
George Moroz <agricolamz@gmail.com>
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/create_subannotation.R
\name{create_subannotation}
\alias{create_subannotation}
\title{Create boundaries in a texgrid tier}
\usage{
create_subannotation(
  textgrid,
  tier = 1,
  new_tier_name = "",
  n_of_annotations = 4,
  each = 1,
  omit_blank = TRUE,
  overwrite = TRUE
)
}
\arguments{
\item{textgrid}{character with a filename or path to the TextGrid}

\item{tier}{value that could be either ordinal number of the tier either name
of the tier}

\item{new_tier_name}{a name of a new created tier}

\item{n_of_annotations}{number of new annotations per annotation to create}

\item{each}{non-negative integer. Each new blank annotation is repeated every
first, second or ... times}

\item{omit_blank}{logical. If TRUE (by dafault) it doesn't create
subannotation for empy annotations.}

\item{overwrite}{logical. If TRUE (by dafault) it overwrites an existing
tier.}
}
\value{
a string that contain TextGrid. If argument write is \code{TRUE},
then no output.
}
\description{
Create boundaries in a texgrid tier
}
\examples{
create_subannotation(system.file("extdata", "test.TextGrid",
  package = "phonfieldwork"
),
tier = 1, overwrite = FALSE
)
}
\author{
George Moroz <agricolamz@gmail.com>
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/remove_textgrid_tier.R
\name{remove_textgrid_tier}
\alias{remove_textgrid_tier}
\title{Remove tier from texgrid}
\usage{
remove_textgrid_tier(textgrid, tier, overwrite = TRUE)
}
\arguments{
\item{textgrid}{character with a filename or path to the TextGrid}

\item{tier}{value that could be either ordinal number of the tier either name
of the tier}

\item{overwrite}{logical. If TRUE (by dafault) it overwrites an existing
tier.}
}
\value{
a string that contain TextGrid. If argument write is \code{TRUE},
then no output.
}
\description{
Remove tier from texgrid
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.r
\name{add_leading_symbols}
\alias{add_leading_symbols}
\title{Create indices padded with zeros}
\usage{
add_leading_symbols(file_names)
}
\arguments{
\item{file_names}{vector of any values.}
}
\value{
A string with numbers padded with leadinng zero.
}
\description{
Create indices padded with zeros. This is important for creating appropriate
for sorting names.
}
\author{
George Moroz <agricolamz@gmail.com>
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/textgrid_to_df.R
\name{textgrid_to_df}
\alias{textgrid_to_df}
\title{TextGrid to dataframe}
\usage{
textgrid_to_df(file_name)
}
\arguments{
\item{file_name}{string with a filename or path to the TextGrid}
}
\value{
a dataframe with columns:  \code{id}, \code{time_start},
\code{time_end} (if it is an interval tier -- the same as the start value),
\code{content}, \code{tier}, \code{tier_name} and \code{source}
}
\description{
Convert Praat TextGrid to a dataframe.
}
\examples{
textgrid_to_df(system.file("extdata", "test.TextGrid",
  package = "phonfieldwork"
))

# this is and example of reading a short .TextGrid format
textgrid_to_df(system.file("extdata", "test_short.TextGrid",
  package = "phonfieldwork"
))
}
\author{
George Moroz <agricolamz@gmail.com>
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/extract_intervals.R
\name{extract_intervals}
\alias{extract_intervals}
\title{Extract intervals}
\usage{
extract_intervals(
  file_name,
  textgrid,
  tier = 1,
  prefix = NULL,
  suffix = NULL,
  autonumber = TRUE,
  path
)
}
\arguments{
\item{file_name}{path to the soundfile}

\item{textgrid}{path to the TextGrid}

\item{tier}{tier number or name that should be used as base for extraction
and names}

\item{prefix}{character vector containing prefix(es) for file names}

\item{suffix}{character vector containing suffix(es) for file names}

\item{autonumber}{if TRUE automatically add number of extracted sound to the
file_name. Prevents from creating a duplicated files and wrong sorting.}

\item{path}{path to the directory where create extracted soundfiles.}
}
\value{
no output
}
\description{
Extract sound according to non-empty annotated intervals from TextGrid and
create soundfiles with correspondent names.
}
\examples{
# create two files in a temprary folder "test_folder"
s <- system.file("extdata", "test.wav", package = "phonfieldwork")
tdir <- tempdir()
file.copy(s, tdir)

# Extract intervals according the TextGrid into the path
extract_intervals(
  file_name = paste0(tdir, "/test.wav"),
  textgrid = system.file("extdata", "test.TextGrid",
    package = "phonfieldwork"
  ),
  path = tdir
)

list.files(tdir)
# [1] "e-2.wav" "s-3.wav" "t-1.wav" "t-4.wav" "test.wav"
}
\author{
George Moroz <agricolamz@gmail.com>
}
