---
title: rfigshare tutorial
author: Carl Boettiger
---


<!--
%\VignetteEngine{knitr::knitr}
%\VignetteIndexEntry{An Introduction to the rfigshare package}
-->

[![Build Status](https://api.travis-ci.org/ropensci/rfigshare.png)](https://travis-ci.org/ropensci/rfigshare)




rfigshare
==========

*An R interface to [FigShare](http://figshare.com)*

* Maintainer: Carl Boettiger, [cboettig](https://github.com/cboettig)
* License: [CC0](http://creativecommons.org/publicdomain/zero/1.0/)
* Contact: Report bugs, questions, or feature requests on the [Issues Tracker](https://github.com/ropensci/rfigshare/issues), or get in touch with us at [info@ropensci.org](mailto:info@ropensci.org)

Installation
------------



```r
require(devtools)
install_github("rfigshare", "ropensci")
```

Getting Started
---------------



# Using rfigshare



```r
require(rfigshare)
```




The first time you use an `rfigshare` function, it will ask you to authenticate online. Just log in and click okay to authenticate rfigshare.  R will allow you to cache your login credentials so that you won't be asked to authenticate again (even between R sessions), as long as you are using the same working directory in future.  

Try a search for an author:



```r
fs_author_search("Boettiger")
```

```
## list()
```



Try creating your own content:



```r
id <- fs_create("Test title", "description of test")
```

```
## Your article has been created! Your id number is 1126334
```


This creates an article with the essential metadata information. We can now upload the dataset to this newly created figshare object using `fs_upload`.  For the purposes of this example, we'll just upload one of R's built-in datasets:



```r
data(mtcars)
write.csv(mtcars, "mtcars.csv")
fs_upload(id, "mtcars.csv")
```


Not that we must pass the id number of our the newly created figshare object as the first argument.  Similar functions let us add additional authors, tags, categories, and links, e.g.



```r
fs_add_tags(id, "demo")
```



Minimal metadata includes title, description, type, and at least one tag and one category.  We can add categories using either the category id or the name of the category, but it must be one of the preset categories available.  We can ask the API for a list of all the categories:



```r
fs_category_list()
```

```
##        id  name                                                      
##   [1,] 1   Biophysics                                                
##   [2,] 4   Biochemistry                                              
##   [3,] 7   Medicine                                                  
##   [4,] 8   Microbiology                                              
##   [5,] 10  Anatomy                                                   
##   [6,] 11  Behavioral Neuroscience                                   
##   [7,] 12  Cell Biology                                              
##   [8,] 13  Genetics                                                  
##   [9,] 14  Molecular Biology                                         
##  [10,] 15  Neuroscience                                              
##  [11,] 16  Physiology                                                
##  [12,] 17  Geography                                                 
##  [13,] 19  Pharmacology                                              
##  [14,] 20  Computer Engineering                                      
##  [15,] 21  Biotechnology                                             
##  [16,] 23  Software Engineering                                      
##  [17,] 24  Evolutionary Biology                                      
##  [18,] 25  Anthropology                                              
##  [19,] 27  Economics                                                 
##  [20,] 28  Paleontology                                              
##  [21,] 29  Geology                                                   
##  [22,] 31  Environmental Chemistry                                   
##  [23,] 32  Geochemistry                                              
##  [24,] 34  Environmental Science                                     
##  [25,] 35  Limnology                                                 
##  [26,] 36  Oceanography                                              
##  [27,] 37  Organic Chemistry                                         
##  [28,] 39  Ecology                                                   
##  [29,] 42  Biological Engineering                                    
##  [30,] 44  Toxicology                                                
##  [31,] 45  Sociology                                                 
##  [32,] 46  Immunology                                                
##  [33,] 47  Stereochemistry                                           
##  [34,] 54  Planetary Geology                                         
##  [35,] 55  Stellar Astronomy                                         
##  [36,] 56  Galactic Astronomy                                        
##  [37,] 57  Cosmology                                                 
##  [38,] 58  Astrophysics                                              
##  [39,] 59  Planetary Science                                         
##  [40,] 61  Developmental Biology                                     
##  [41,] 62  Marine Biology                                            
##  [42,] 63  Parasitology                                              
##  [43,] 64  Cancer                                                    
##  [44,] 65  Botany                                                    
##  [45,] 66  Crystallography                                           
##  [46,] 69  Inorganic Chemistry                                       
##  [47,] 70  Molecular Physics                                         
##  [48,] 71  Nuclear Chemistry                                         
##  [49,] 73  Radiochemistry                                            
##  [50,] 75  Supramolecular Chemistry                                  
##  [51,] 76  Theoretical Computer Science                              
##  [52,] 77  Applied Computer Science                                  
##  [53,] 78  Atmospheric Sciences                                      
##  [54,] 79  Geophysics                                                
##  [55,] 80  Hydrology                                                 
##  [56,] 82  Mineralogy                                                
##  [57,] 84  Paleoclimatology                                          
##  [58,] 85  Palynology                                                
##  [59,] 86  Physical Geography                                        
##  [60,] 87  Soil Science                                              
##  [61,] 88  Agricultural Engineering                                  
##  [62,] 89  Aerospace Engineering                                     
##  [63,] 90  Genetic Engineering                                       
##  [64,] 91  Mechanical Engineering                                    
##  [65,] 92  Nuclear Engineering                                       
##  [66,] 94  Art                                                       
##  [67,] 95  Design                                                    
##  [68,] 97  Law                                                       
##  [69,] 98  Literature                                                
##  [70,] 99  Performing Arts                                           
##  [71,] 100 Philosophy                                                
##  [72,] 102 Algebra                                                   
##  [73,] 103 Geometry                                                  
##  [74,] 104 Probability                                               
##  [75,] 105 Statistics                                                
##  [76,] 106 Science Policy                                            
##  [77,] 107 Applied Physics                                           
##  [78,] 108 Atomic Physics                                            
##  [79,] 109 Computational Physics                                     
##  [80,] 110 Condensed Matter Physics                                  
##  [81,] 111 Mechanics                                                 
##  [82,] 112 Particle Physics                                          
##  [83,] 113 Plasma Physics                                            
##  [84,] 114 Quantum Mechanics                                         
##  [85,] 115 Solid Mechanics                                           
##  [86,] 116 Thermodynamics                                            
##  [87,] 117 Entropy                                                   
##  [88,] 118 General Relativity                                        
##  [89,] 119 M-Theory                                                  
##  [90,] 120 Special Relativity                                        
##  [91,] 122 Mental Health                                             
##  [92,] 125 Bioinformatics                                            
##  [93,] 126 History                                                   
##  [94,] 127 Archaeology                                               
##  [95,] 128 Hematology                                                
##  [96,] 129 Education                                                 
##  [97,] 130 Survey Results                                            
##  [98,] 131 Anesthesiology                                            
##  [99,] 132 Infectious Diseases                                       
## [100,] 133 Plant Biology                                             
## [101,] 134 Virology                                                  
## [102,] 135 Computational  Biology                                    
## [103,] 136 NMR Spectroscopy                                          
## [104,] 137 Cheminformatics                                           
## [105,] 138 Numerical Analysis                                        
## [106,] 139 Pathology                                                 
## [107,] 140 Cardiology                                                
## [108,] 141 Computational Chemistry                                   
## [109,] 143 Solid Earth Sciences                                      
## [110,] 144 Climate Science                                           
## [111,] 145 Solar System, Solar Physics, Planets and Exoplanets       
## [112,] 146 Space Science                                             
## [113,] 147 Stars, Variable Stars                                     
## [114,] 148 Instrumentation, Techniques, and Astronomical Observations
## [115,] 149 Interstellar and Intergalactic Matter                     
## [116,] 150 Extragalactic Astronomy                                   
## [117,] 151 Biomarkers                                                
## [118,] 152 Pathogenesis                                              
## [119,] 153 Health Care                                               
## [120,] 154 Diseases                                                  
## [121,] 155 Stem Cells                                                
## [122,] 156 Systems Biology                                           
## [123,] 157 Structural Biology                                        
## [124,] 158 Biological Techniques                                     
## [125,] 159 Zoology                                                   
## [126,] 160 Digital Humanities                                        
## [127,] 161 Disability Studies                                        
## [128,] 162 Drama                                                     
## [129,] 163 Entertainment                                             
## [130,] 164 Environmental Humanities                                  
## [131,] 165 Ethnic Studies                                            
## [132,] 166 Gender studies                                            
## [133,] 167 Language                                                  
## [134,] 168 Linguistics                                               
## [135,] 169 Media Studies                                             
## [136,] 170 Museology                                                 
## [137,] 171 Religious Studies                                         
## [138,] 172 Rhetoric                                                  
## [139,] 173 Applied Psychology                                        
## [140,] 174 Clinical Psychology                                       
## [141,] 175 Developmental and Educational Psychology                  
## [142,] 176 Neuroscience and Physiological Psychology                 
## [143,] 177 Organizational Behavioral Psychology                      
## [144,] 178 Personality, Social and Criminal Psychology               
## [145,] 179 Artificial Intelligence and Image Processing              
## [146,] 180 Computation Theory and Mathematics                        
## [147,] 181 Computer Software                                         
## [148,] 182 Data Format                                               
## [149,] 183 Distributed Computing                                     
## [150,] 184 Information Systems                                       
## [151,] 185 Library and Information Studies
```


And we can add the category or categories we like,



```r
fs_add_categories(id, c("Education", "Software Engineering"))
```



The file we have created remains saved as a draft until we publish it, either publicly or privately.  Note that once a file is declared public, it cannot be made private or deleted.  Let's release this dataset privately:



```r
fs_make_private(id)
```

```
## Response [http://api.figshare.com/v1/my_data/articles/1126334/action/make_private]
##   Status: 200
##   Content-type: application/json; charset=UTF-8
##   Size: 48 B
## {"success": "Article status changed to Private"}
```


We can now share the dataset with collaborators by way of the private url.  

### The quick and easy way

The `rfigshare` package will let you create a new figshare article with additional authors, tags, categories, etc in a single command usnig the `fs_new_article` function.  The essential metadata `title`, `description` and `type` are required, but any other information we omit at this stage can be added later.  If we set `visibility` to private or public, the article is published on figshare immediately.  



```r
data(mtcars)
write.csv(mtcars,"mtcars.csv")
id <- fs_new_article(title="A Test of rfigshare", 
                     description="This is a test", 
                     type="dataset", 
                     authors=c("Karthik Ram", "Scott Chamberlain"), 
                     tags=c("ecology", "openscience"), 
                     categories="Ecology", 
                     links="http://ropensci.org", 
                     files="mtcars.csv",
                     visibility="private")
```

```
## Your article has been created! Your id number is 1126335
```

```r
unlink("mtcars.csv") # clean up
```


# Examining Data on Figshare

We can view all available metadata of a figshare object. 



```r
fs_details(id)
```

```
## article_id: 1.1263e+06
## title: A Test of rfigshare
## master_publisher_id: 0.0e+00
## defined_type: dataset
## status: Private
## version: 1.0
## published_date: 06:25, Aug 03, 2014
## description: This is a test
## description_nohtml: This is a test
## total_size: 1.70 KB
## authors:
## - first_name: Ropensci
##   last_name: Testaccount
##   id: 4.3142e+05
##   full_name: Ropensci TestAccount
## - first_name: Karthik
##   last_name: Ram
##   id: 9.7306e+04
##   full_name: Karthik Ram
## - first_name: Scott
##   last_name: Chamberlain
##   id: 9.6554e+04
##   full_name: Scott Chamberlain
## tags:
## - id: 1.9539e+05
##   name: Published using rfigshare
## - id: 4.7864e+04
##   name: openscience
## - id: 1.1917e+04
##   name: ecology
## categories:
## - id: 39.0
##   name: Ecology
## files:
## - size: 2 KB
##   thumb: ~
##   id: 1.6204e+06
##   mime_type: text/plain
##   name: mtcars.csv
## links:
## - link: http://ropensci.org
##   id: 673.0
```

You can see all of the files you have (Currently only up to 10):



```r
mine <- fs_browse()
mine[1:2]
```

```
## [[1]]
## article_id: 1.1263e+06
## title: A Test of rfigshare
## master_publisher_id: 0.0e+00
## defined_type: dataset
## status: Private
## version: 1.0
## published_date: 06:25, Aug 03, 2014
## description: This is a test
## description_nohtml: This is a test
## total_size: 1.70 KB
## authors:
## - first_name: Ropensci
##   last_name: Testaccount
##   id: 4.3142e+05
##   full_name: Ropensci TestAccount
## - first_name: Karthik
##   last_name: Ram
##   id: 9.7306e+04
##   full_name: Karthik Ram
## - first_name: Scott
##   last_name: Chamberlain
##   id: 9.6554e+04
##   full_name: Scott Chamberlain
## tags:
## - id: 1.9539e+05
##   name: Published using rfigshare
## - id: 4.7864e+04
##   name: openscience
## - id: 1.1917e+04
##   name: ecology
## categories:
## - id: 39.0
##   name: Ecology
## files:
## - size: 2 KB
##   thumb: ~
##   id: 1.6204e+06
##   mime_type: text/plain
##   name: mtcars.csv
## links:
## - link: http://ropensci.org
##   id: 673.0
## 
## [[2]]
## article_id: 1.1263e+06
## title: Test title
## master_publisher_id: 0.0e+00
## defined_type: dataset
## status: Private
## version: 1.0
## published_date: 06:25, Aug 03, 2014
## description: description of test
## description_nohtml: description of test
## total_size: 1.70 KB
## authors:
## - first_name: Ropensci
##   last_name: Testaccount
##   id: 4.3142e+05
##   full_name: Ropensci TestAccount
## tags:
## - id: 1.5681e+04
##   name: demo
## - id: 1.9539e+05
##   name: Published using rfigshare
## categories:
## - id: 23.0
##   name: Software Engineering
## - id: 129.0
##   name: Education
## files:
## - size: 2 KB
##   thumb: ~
##   id: 1.6204e+06
##   mime_type: text/plain
##   name: mtcars.csv
## links: []
```

Note that we can easily grab the ids with the wrapper function `fs_ids`:




```r
fs_ids(mine)
```

```
##  [1] 1126335 1126334 1126332 1126329 1126328 1126324 1126322 1126321
##  [9] 1126318 1126317
```


We can delete unwanted files that are not public with `fs_delete`:  



```r
fs_delete(id)
```

To cite package `rfigshare` in publications use:



```r
citation("rfigshare")
```

```
## 
## To cite package 'rfigshare' in publications use:
## 
##   Carl Boettiger, Scott Chamberlain, Karthik Ram and Edmund Hart
##   (2014). rfigshare: an R interface to figshare.com.. R package
##   version 0.3-1. http://CRAN.R-project.org/package=rfigshare
## 
## A BibTeX entry for LaTeX users is
## 
##   @Manual{,
##     title = {rfigshare: an R interface to figshare.com.},
##     author = {Carl Boettiger and Scott Chamberlain and Karthik Ram and Edmund Hart},
##     year = {2014},
##     note = {R package version 0.3-1},
##     url = {http://CRAN.R-project.org/package=rfigshare},
##   }
```


[![ropensci footer](http://ropensci.org/public_images/github_footer.png)](http://ropensci.org)
Dear CRAN maintainers,

This release fixes the title format and should address a compatibility issue with the httr package dependency prior to the release of httr 1.0.
I have now also addressed the issues in DESCRIPTION file, my apologies for not addressing these in the original update.

Sincerely,

Carl Boettiger
rOpenSci, Figshare and knitr.
========================================================

This is a tutorial of how you can create reproducible documents using [knitr]( http://yihui.name/knitr/) and pandoc, and seamlessly upload them to [figshare](http;//www.figshare.com) attaching a citable DOI to them. This document will walk you through the process of creating a document in knitr and uploading a compiled PDF to Figshare.  I make the following assumptions about your knowledge:

* You have set-up an account at [Figshare.com](http://www.figshare.com)

* You have installed [rfigshare](http://cran.r-project.org/web/packages/rfigshare/index.html), the [knitr](http://yihui.name/knitr/) package and are familiar with the concepts of knitr or sweave, as well as [pandoc](http://johnmacfarlane.net/pandoc/), for conversion to pdf (although this is only necessary if you want to convert your document to a pdf)

* You have successfully set-up the your credentials for rfigshare.  If not go to our [tutorial](http://github.com/ropensci/rfigshare/blob/master/inst/doc/tutorial.md) and make sure your credentials are properly set.


The goal of this document is to demonstrate how one could carry out a project using tools from [rOpenSci](http://www.ropensci.org/), knitr, and share the results on figshare in one continuous workflow.  To do this I'll be using a tutorial from one of our packages, [treebase](http://cran.r-project.org/web/packages/treebase/), which allows you to download trees from [TreebaseWEB](http://treebase.org/treebase-web/home.html;jsessionid=A258F89FBF584F44E0CDB740B8ECF3A8)


First I'll turn the cache on.

```r
opts_chunk$set(cache = TRUE, autodep = TRUE)
dep_auto()
```

Then you'll want to download some data, and maybe make a plot, and say some things about how great your plot is.

```r
library(treebase)
tree <- search_treebase("Derryberry", "author")[[1]]
# plotting only part of the tree because it's so large
plot.phylo(tree, y.lim = c(0, 20))
```

![My amazing tree!](figure/unnamed-chunk-1.png) 


Once you've made all your plots, and said all you want to say it's time to convert your document, and then create a new article using `fs_new_article()`



```r
library(knitr)
library(rfigshare)
options(FigshareKey = "XXXXXXXX")
options(FigsharePrivateKey = "XXXXXXXX")
options(FigshareToken = "XXXXXXXX")
options(FigsharePrivateToken = "XXXXXXXX")

#knit document to pandoc markdown
knit("rfigtutorial.Rmd")
#convert to pdf
system("pandoc -S rfigtutorial.md -o rfigtutorial.pdf")

id <- fs_new_article(title="An rfigshare tutorial", 
                      description="How to create a document in knitr and 
                      upload it to figshare.com", 
                     
                      type="paper", 
                      authors=c("Edmund Hart"), 
                      tags=c("ecology", "openscience"), 
                      categories="Ecology", 
                      links="http://emhart.github.com", 
                      files="rfigtutorial.pdf",
                      visibility="draft")
```


The main advantage of this approach is that manuscripts can be worked on from within the R environment and then seemlessly uploaded to figshare. Also it's best practice to store your key values in your `.Rprofile` so I would reccomend file Be sure to run `fs_make_public(id)` when you're ready to make your article public.  


---
title: rfigshare - Share and explore figures, data, and publications on FigShare using R

---

![](http://farm9.staticflickr.com/8180/7950489358_ea902bdaae_o.png)







# Obtaining your API keys

There is a nice video introduction to creating applications for the API on the [figshare blog](http://figshare.com/blog/figshare_API_available_to_all/48).  The following tutorial provides a simple walkthrough of how to go about getting your figshare API keys set up so that you can use the `rfigshare` package.  


Create a user account on [FigShare](http://figshare.com) and log in.  From your homepage, select "Applications" from the drop-down menu,

![](http://farm9.staticflickr.com/8171/7950489558_5172515057_o.png)

Create a new application:

![](http://farm9.staticflickr.com/8038/7950490158_7feaf680bd_o.png)


Enter in the following information: 

![](http://farm9.staticflickr.com/8305/7950490562_02846cea92_o.png)

Then navigate over to the permissions tab.  To get the most out of `rfigshare` you'll want to enable all permissions:

![](http://farm9.staticflickr.com/8448/7950491064_c3820e62bd_o.png)

Save the new settings, and then open the application again (View/Edit menu) and click on the "Access Codes" tab.

![](http://farm9.staticflickr.com/8308/7950491470_621da9c5d1_o.png)

Record each if the keys into R as follows.  You might want to put this bit of R code into your `.Rprofile` to avoid entering it each time in the future:

```r
options(FigshareKey = "qMDabXXXXXXXXXXXXXXXXX")
options(FigsharePrivateKey = "zQXXXXXXXXXXXXXXXXXXXX")
options(FigshareToken = "SrpxabQXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX")
options(FigsharePrivateToken = "yqXXXXXXXXXXXXXXXXXXXX")
```


## Installing the R package

Now that we have the API credentials in place, actually using `rfigshare` is pretty easy.  Install the latest version of package directly from Github using: 

```r
require(devtools)
install_github("rfigshare", "ropensci")
```

# Using rfigshare


Try authenticating with your credentials:


```r
require(rfigshare)
fs_auth()
```



Try a search for an author:


```r
fs_author_search("Boettiger")
```

```
text_content() deprecated. Use parsed_content(x, as = 'parsed')
```

```
[[1]]
[[1]]$id
[1] "96387"

[[1]]$fname
[1] "Carl"

[[1]]$lname
[1] "Boettiger"

[[1]]$full_name
[1] "Carl Boettiger"

[[1]]$job_title
[1] ""

[[1]]$description
[1] ""

[[1]]$facebook
[1] ""

[[1]]$twitter
[1] ""

[[1]]$active
[1] 1

```


Try creating your own content:


```r
id <- fs_create("Test title", "description of test", "dataset")
```

```
text_content() deprecated. Use parsed_content(x, as = 'parsed')
```

```
Your article has been created! Your id number is 105137
```


This creates an article with the essential metadata information. We can now upload the dataset to this newly created figshare object using `fs_upload`.  For the purposes of this example, we'll just upload one of R's built-in datasets:


```r
data(mtcars)
write.csv(mtcars, "mtcars.csv")

fs_upload(id, "mtcars.csv")
```


Not that we must pass the id number of our the newly created figshare object as the first argument.  Similar functions let us add additional authors, tags, categories, and links, e.g.


```r
fs_add_tags(id, "demo")
```



Minimal metadata includes title, description, type, and at least one tag and one category.  We can add categories using either the category id or the name of the category, but it must be one of the preset categories available.  We can ask the API for a list of all the categories:


```r
fs_category_list()
```


<!-- html table generated in R 2.15.2 by xtable 1.7-0 package -->
<!-- Wed Dec 19 22:02:19 2012 -->
<TABLE border=1>
<TR> <TH>  </TH> <TH> id </TH> <TH> name </TH>  </TR>
  <TR> <TD align="right"> 1 </TD> <TD> 89 </TD> <TD> Aerospace Engineering </TD> </TR>
  <TR> <TD align="right"> 2 </TD> <TD> 88 </TD> <TD> Agricultural Engineering </TD> </TR>
  <TR> <TD align="right"> 3 </TD> <TD> 102 </TD> <TD> Algebra </TD> </TR>
  <TR> <TD align="right"> 4 </TD> <TD> 10 </TD> <TD> Anatomy </TD> </TR>
  <TR> <TD align="right"> 5 </TD> <TD> 25 </TD> <TD> Anthropology </TD> </TR>
  <TR> <TD align="right"> 6 </TD> <TD> 77 </TD> <TD> Applied Computer Science </TD> </TR>
  <TR> <TD align="right"> 7 </TD> <TD> 107 </TD> <TD> Applied Physics </TD> </TR>
  <TR> <TD align="right"> 8 </TD> <TD> 127 </TD> <TD> Archaeology </TD> </TR>
  <TR> <TD align="right"> 9 </TD> <TD> 94 </TD> <TD> Art </TD> </TR>
  <TR> <TD align="right"> 10 </TD> <TD> 58 </TD> <TD> Astrophysics </TD> </TR>
  <TR> <TD align="right"> 11 </TD> <TD> 78 </TD> <TD> Atmospheric Sciences </TD> </TR>
  <TR> <TD align="right"> 12 </TD> <TD> 108 </TD> <TD> Atomic Physics </TD> </TR>
  <TR> <TD align="right"> 13 </TD> <TD> 11 </TD> <TD> Behavioral neuroscience </TD> </TR>
  <TR> <TD align="right"> 14 </TD> <TD> 4 </TD> <TD> Biochemistry </TD> </TR>
  <TR> <TD align="right"> 15 </TD> <TD> 125 </TD> <TD> Bioinformatics </TD> </TR>
  <TR> <TD align="right"> 16 </TD> <TD> 42 </TD> <TD> Biological engineering </TD> </TR>
  <TR> <TD align="right"> 17 </TD> <TD> 1 </TD> <TD> Biophysics </TD> </TR>
  <TR> <TD align="right"> 18 </TD> <TD> 21 </TD> <TD> Biotechnology </TD> </TR>
  <TR> <TD align="right"> 19 </TD> <TD> 65 </TD> <TD> Botany </TD> </TR>
  <TR> <TD align="right"> 20 </TD> <TD> 64 </TD> <TD> Cancer </TD> </TR>
  <TR> <TD align="right"> 21 </TD> <TD> 12 </TD> <TD> Cell biology </TD> </TR>
  <TR> <TD align="right"> 22 </TD> <TD> 109 </TD> <TD> Computational Physics </TD> </TR>
  <TR> <TD align="right"> 23 </TD> <TD> 20 </TD> <TD> Computer Engineering </TD> </TR>
  <TR> <TD align="right"> 24 </TD> <TD> 110 </TD> <TD> Condensed Matter Physics </TD> </TR>
  <TR> <TD align="right"> 25 </TD> <TD> 57 </TD> <TD> Cosmology </TD> </TR>
  <TR> <TD align="right"> 26 </TD> <TD> 66 </TD> <TD> Crystallography </TD> </TR>
  <TR> <TD align="right"> 27 </TD> <TD> 95 </TD> <TD> Design </TD> </TR>
  <TR> <TD align="right"> 28 </TD> <TD> 61 </TD> <TD> Developmental Biology </TD> </TR>
  <TR> <TD align="right"> 29 </TD> <TD> 39 </TD> <TD> Ecology </TD> </TR>
  <TR> <TD align="right"> 30 </TD> <TD> 27 </TD> <TD> Economics </TD> </TR>
  <TR> <TD align="right"> 31 </TD> <TD> 129 </TD> <TD> Education </TD> </TR>
  <TR> <TD align="right"> 32 </TD> <TD> 117 </TD> <TD> Entropy </TD> </TR>
  <TR> <TD align="right"> 33 </TD> <TD> 31 </TD> <TD> Environmental chemistry </TD> </TR>
  <TR> <TD align="right"> 34 </TD> <TD> 34 </TD> <TD> Environmental science </TD> </TR>
  <TR> <TD align="right"> 35 </TD> <TD> 24 </TD> <TD> Evolutionary biology </TD> </TR>
  <TR> <TD align="right"> 36 </TD> <TD> 56 </TD> <TD> Galactic Astronomy </TD> </TR>
  <TR> <TD align="right"> 37 </TD> <TD> 118 </TD> <TD> General Relativity </TD> </TR>
  <TR> <TD align="right"> 38 </TD> <TD> 90 </TD> <TD> Genetic Engineering </TD> </TR>
  <TR> <TD align="right"> 39 </TD> <TD> 13 </TD> <TD> Genetics </TD> </TR>
  <TR> <TD align="right"> 40 </TD> <TD> 32 </TD> <TD> Geochemistry </TD> </TR>
  <TR> <TD align="right"> 41 </TD> <TD> 17 </TD> <TD> Geography </TD> </TR>
  <TR> <TD align="right"> 42 </TD> <TD> 29 </TD> <TD> Geology </TD> </TR>
  <TR> <TD align="right"> 43 </TD> <TD> 103 </TD> <TD> Geometry </TD> </TR>
  <TR> <TD align="right"> 44 </TD> <TD> 79 </TD> <TD> Geophysics </TD> </TR>
  <TR> <TD align="right"> 45 </TD> <TD> 128 </TD> <TD> Hematology </TD> </TR>
  <TR> <TD align="right"> 46 </TD> <TD> 126 </TD> <TD> History </TD> </TR>
  <TR> <TD align="right"> 47 </TD> <TD> 80 </TD> <TD> Hydrology </TD> </TR>
  <TR> <TD align="right"> 48 </TD> <TD> 46 </TD> <TD> Immunology </TD> </TR>
  <TR> <TD align="right"> 49 </TD> <TD> 69 </TD> <TD> Inorganic Chemistry </TD> </TR>
  <TR> <TD align="right"> 50 </TD> <TD> 97 </TD> <TD> Law </TD> </TR>
  <TR> <TD align="right"> 51 </TD> <TD> 35 </TD> <TD> Limnology </TD> </TR>
  <TR> <TD align="right"> 52 </TD> <TD> 98 </TD> <TD> Literature </TD> </TR>
  <TR> <TD align="right"> 53 </TD> <TD> 119 </TD> <TD> M-Theory </TD> </TR>
  <TR> <TD align="right"> 54 </TD> <TD> 62 </TD> <TD> Marine Biology </TD> </TR>
  <TR> <TD align="right"> 55 </TD> <TD> 91 </TD> <TD> Mechanical Engineering </TD> </TR>
  <TR> <TD align="right"> 56 </TD> <TD> 111 </TD> <TD> Mechanics </TD> </TR>
  <TR> <TD align="right"> 57 </TD> <TD> 7 </TD> <TD> Medicine </TD> </TR>
  <TR> <TD align="right"> 58 </TD> <TD> 122 </TD> <TD> Mental Health </TD> </TR>
  <TR> <TD align="right"> 59 </TD> <TD> 8 </TD> <TD> Microbiology </TD> </TR>
  <TR> <TD align="right"> 60 </TD> <TD> 82 </TD> <TD> Mineralogy </TD> </TR>
  <TR> <TD align="right"> 61 </TD> <TD> 14 </TD> <TD> Molecular biology </TD> </TR>
  <TR> <TD align="right"> 62 </TD> <TD> 70 </TD> <TD> Molecular Physics </TD> </TR>
  <TR> <TD align="right"> 63 </TD> <TD> 15 </TD> <TD> Neuroscience </TD> </TR>
  <TR> <TD align="right"> 64 </TD> <TD> 71 </TD> <TD> Nuclear Chemistry </TD> </TR>
  <TR> <TD align="right"> 65 </TD> <TD> 92 </TD> <TD> Nuclear Engineering </TD> </TR>
  <TR> <TD align="right"> 66 </TD> <TD> 36 </TD> <TD> Oceanography </TD> </TR>
  <TR> <TD align="right"> 67 </TD> <TD> 37 </TD> <TD> Organic chemistry </TD> </TR>
  <TR> <TD align="right"> 68 </TD> <TD> 84 </TD> <TD> Paleoclimatology </TD> </TR>
  <TR> <TD align="right"> 69 </TD> <TD> 28 </TD> <TD> Paleontology </TD> </TR>
  <TR> <TD align="right"> 70 </TD> <TD> 85 </TD> <TD> Palynology </TD> </TR>
  <TR> <TD align="right"> 71 </TD> <TD> 63 </TD> <TD> Parasitology </TD> </TR>
  <TR> <TD align="right"> 72 </TD> <TD> 112 </TD> <TD> Particle Physics </TD> </TR>
  <TR> <TD align="right"> 73 </TD> <TD> 99 </TD> <TD> Performing Arts </TD> </TR>
  <TR> <TD align="right"> 74 </TD> <TD> 19 </TD> <TD> Pharmacology </TD> </TR>
  <TR> <TD align="right"> 75 </TD> <TD> 100 </TD> <TD> Philosophy </TD> </TR>
  <TR> <TD align="right"> 76 </TD> <TD> 86 </TD> <TD> Physical Geography </TD> </TR>
  <TR> <TD align="right"> 77 </TD> <TD> 16 </TD> <TD> Physiology </TD> </TR>
  <TR> <TD align="right"> 78 </TD> <TD> 54 </TD> <TD> Planetary Geology </TD> </TR>
  <TR> <TD align="right"> 79 </TD> <TD> 59 </TD> <TD> Planetary Science </TD> </TR>
  <TR> <TD align="right"> 80 </TD> <TD> 113 </TD> <TD> Plasma Physics </TD> </TR>
  <TR> <TD align="right"> 81 </TD> <TD> 104 </TD> <TD> Probability </TD> </TR>
  <TR> <TD align="right"> 82 </TD> <TD> 114 </TD> <TD> Quantum Mechanics </TD> </TR>
  <TR> <TD align="right"> 83 </TD> <TD> 73 </TD> <TD> Radiochemistry </TD> </TR>
  <TR> <TD align="right"> 84 </TD> <TD> 106 </TD> <TD> Science Policy </TD> </TR>
  <TR> <TD align="right"> 85 </TD> <TD> 45 </TD> <TD> Sociology </TD> </TR>
  <TR> <TD align="right"> 86 </TD> <TD> 23 </TD> <TD> Software Engineering </TD> </TR>
  <TR> <TD align="right"> 87 </TD> <TD> 87 </TD> <TD> Soil Science </TD> </TR>
  <TR> <TD align="right"> 88 </TD> <TD> 115 </TD> <TD> Solid Mechanics </TD> </TR>
  <TR> <TD align="right"> 89 </TD> <TD> 120 </TD> <TD> Special Relativity </TD> </TR>
  <TR> <TD align="right"> 90 </TD> <TD> 105 </TD> <TD> Statistics </TD> </TR>
  <TR> <TD align="right"> 91 </TD> <TD> 55 </TD> <TD> Stellar Astronomy </TD> </TR>
  <TR> <TD align="right"> 92 </TD> <TD> 47 </TD> <TD> Stereochemistry </TD> </TR>
  <TR> <TD align="right"> 93 </TD> <TD> 75 </TD> <TD> Supramolecular Chemistry </TD> </TR>
  <TR> <TD align="right"> 94 </TD> <TD> 76 </TD> <TD> Theoretical Computer Science </TD> </TR>
  <TR> <TD align="right"> 95 </TD> <TD> 116 </TD> <TD> Thermodynamics </TD> </TR>
  <TR> <TD align="right"> 96 </TD> <TD> 44 </TD> <TD> Toxicology </TD> </TR>
   </TABLE>


And we can add the category or categories we like,


```r
fs_add_categories(id, c("Education", "Software Engineering"))
```



The file we have created remains saved as a draft until we publish it, either publicly or privately.  Note that once a file is declared public, it cannot be made private or deleted.  Let's release this dataset privately:


```r
fs_make_private(id)
```

```
Response [http://api.figshare.com/v1/my_data/articles/105137/action/make_private]
  Status: 200
  Content-type: application/json; charset=UTF-8
{"success": "Article status changed to Private"} 
```


We can now share the dataset with collaborators by way of the private url.  

### The quick and easy way

The `rfigshare` package will let you create a new figshare article with additional authors, tags, categories, etc in a single command usnig the `fs_new_article` function.  The essential metadata `title`, `description` and `type` are required, but any other information we omit at this stage can be added later.  If we set `visibility` to private or public, the article is published on figshare immediately.  


```r
id <- fs_new_article(title="A Test of rfigshare", 
                     description="This is a test", 
                     type="figure", 
                     authors=c("Karthik Ram", "Scott Chamberlain"), 
                     tags=c("ecology", "openscience"), 
                     categories="Ecology", 
                     links="http://ropensci.org", 
                     files="figure/rfigshare.png",
                     visibility="private")
```

```
text_content() deprecated. Use parsed_content(x, as = 'parsed')
```

```
Your article has been created! Your id number is 105138
```

```
text_content() deprecated. Use parsed_content(x, as = 'parsed')
```

```
text_content() deprecated. Use parsed_content(x, as = 'parsed')
```

```
found ids for all authors
```


# Examining Data on Figshare

We can view all available metadata of a figshare object.  If the object is not public (e.g. draft or private), we have to add the `mine=TRUE` option


```r
fs_details(id, mine=TRUE)
```

```
text_content() deprecated. Use parsed_content(x, as = 'parsed')
```

```
Article ID : 105138
Article type : figure
DOI : 
Title : A Test of rfigshare
Description : This is a test
Shares : 
Views : 
Downloads : 
Owner : 
Authors : Carl Boettiger, Karthik Ram, Scott Chamberlain
Tags : openscience, ecology
Categories : Ecology
File names : rfigshare.png
Links : http://ropensci.org
```


You can see all of the files you have:


```r
fs_browse(mine=TRUE)
```

```
text_content() deprecated. Use parsed_content(x, as = 'parsed')
```

```
[[1]]
[[1]]$article_id
[1] 105138

[[1]]$title
[1] "A Test of rfigshare"

[[1]]$defined_type
[1] "figure"

[[1]]$status
[1] "Private"

[[1]]$version
[1] 1

[[1]]$published_date
[1] "05:18, Dec 20, 2012"

[[1]]$description
[1] "This is a test"

[[1]]$description_nohtml
[1] "This is a test"

[[1]]$total_size
[1] "17.78 KB"

[[1]]$authors
[[1]]$authors[[1]]
[[1]]$authors[[1]]$first_name
[1] "Carl"

[[1]]$authors[[1]]$last_name
[1] "Boettiger"

[[1]]$authors[[1]]$id
[1] 96387

[[1]]$authors[[1]]$full_name
[1] "Carl Boettiger"


[[1]]$authors[[2]]
[[1]]$authors[[2]]$first_name
[1] "Karthik"

[[1]]$authors[[2]]$last_name
[1] "Ram"

[[1]]$authors[[2]]$id
[1] 97306

[[1]]$authors[[2]]$full_name
[1] "Karthik Ram"


[[1]]$authors[[3]]
[[1]]$authors[[3]]$first_name
[1] "Scott"

[[1]]$authors[[3]]$last_name
[1] "Chamberlain"

[[1]]$authors[[3]]$id
[1] 96554

[[1]]$authors[[3]]$full_name
[1] "Scott Chamberlain"



[[1]]$tags
[[1]]$tags[[1]]
[[1]]$tags[[1]]$id
[1] 47864

[[1]]$tags[[1]]$name
[1] "openscience"


[[1]]$tags[[2]]
[[1]]$tags[[2]]$id
[1] 11917

[[1]]$tags[[2]]$name
[1] "ecology"



[[1]]$categories
[[1]]$categories[[1]]
[[1]]$categories[[1]]$id
[1] 39

[[1]]$categories[[1]]$name
[1] "Ecology"



[[1]]$files
[[1]]$files[[1]]
[[1]]$files[[1]]$size
[1] "18 KB"

[[1]]$files[[1]]$id
[1] 233421

[[1]]$files[[1]]$mime_type
[1] "image/png"

[[1]]$files[[1]]$name
[1] "rfigshare.png"



[[1]]$links
[[1]]$links[[1]]
[[1]]$links[[1]]$link
[1] "http://ropensci.org"

[[1]]$links[[1]]$id
[1] 673




[[2]]
[[2]]$article_id
[1] 105137

[[2]]$title
[1] "Test title"

[[2]]$defined_type
[1] "dataset"

[[2]]$status
[1] "Private"

[[2]]$version
[1] 1

[[2]]$published_date
[1] "05:18, Dec 20, 2012"

[[2]]$description
[1] "description of test"

[[2]]$description_nohtml
[1] "description of test"

[[2]]$total_size
[1] "1.70 KB"

[[2]]$authors
[[2]]$authors[[1]]
[[2]]$authors[[1]]$first_name
[1] "Carl"

[[2]]$authors[[1]]$last_name
[1] "Boettiger"

[[2]]$authors[[1]]$id
[1] 96387

[[2]]$authors[[1]]$full_name
[1] "Carl Boettiger"



[[2]]$tags
[[2]]$tags[[1]]
[[2]]$tags[[1]]$id
[1] 15681

[[2]]$tags[[1]]$name
[1] "demo"



[[2]]$categories
[[2]]$categories[[1]]
[[2]]$categories[[1]]$id
[1] 23

[[2]]$categories[[1]]$name
[1] "Software Engineering"


[[2]]$categories[[2]]
[[2]]$categories[[2]]$id
[1] 129

[[2]]$categories[[2]]$name
[1] "Education"



[[2]]$files
[[2]]$files[[1]]
[[2]]$files[[1]]$size
[1] "2 KB"

[[2]]$files[[1]]$id
[1] 233420

[[2]]$files[[1]]$mime_type
[1] "text/plain"

[[2]]$files[[1]]$name
[1] "mtcars.csv"



[[2]]$links
list()


[[3]]
[[3]]$article_id
[1] 97653

[[3]]$title
[1] "RFigshare Tutorial"

[[3]]$defined_type
[1] "paper"

[[3]]$status
[1] "Public"

[[3]]$version
[1] 2

[[3]]$published_date
[1] "17:16, Nov 16, 2012"

[[3]]$description
[1] "<p>A tuturial on how to setup and post material to figshare from R</p>"

[[3]]$description_nohtml
[1] "A tuturial on how to setup and post material to figshare from R"

[[3]]$total_size
[1] "16.21 KB"

[[3]]$authors
[[3]]$authors[[1]]
[[3]]$authors[[1]]$first_name
[1] "Edmund"

[[3]]$authors[[1]]$last_name
[1] "Hart"

[[3]]$authors[[1]]$id
[1] 98137

[[3]]$authors[[1]]$full_name
[1] "Edmund Hart"


[[3]]$authors[[2]]
[[3]]$authors[[2]]$first_name
[1] "Carl"

[[3]]$authors[[2]]$last_name
[1] "Boettiger"

[[3]]$authors[[2]]$id
[1] 96387

[[3]]$authors[[2]]$full_name
[1] "Carl Boettiger"



[[3]]$tags
[[3]]$tags[[1]]
[[3]]$tags[[1]]$id
[1] 47864

[[3]]$tags[[1]]$name
[1] "openscience"


[[3]]$tags[[2]]
[[3]]$tags[[2]]$id
[1] 11917

[[3]]$tags[[2]]$name
[1] "ecology"



[[3]]$categories
[[3]]$categories[[1]]
[[3]]$categories[[1]]$id
[1] 39

[[3]]$categories[[1]]$name
[1] "Ecology"



[[3]]$files
[[3]]$files[[1]]
[[3]]$files[[1]]$download_url
[1] "http://files.figshare.com/223958/tutorial.md"

[[3]]$files[[1]]$size
[1] "17 KB"

[[3]]$files[[1]]$id
[1] 223958

[[3]]$files[[1]]$mime_type
[1] "text/plain"

[[3]]$files[[1]]$name
[1] "tutorial.md"



[[3]]$links
[[3]]$links[[1]]
[[3]]$links[[1]]$link
[1] "http://ropensci.org"

[[3]]$links[[1]]$id
[1] 673




[[4]]
[[4]]$article_id
[1] 97500

[[4]]$title
[1] "Exit Seminar to Center for Population Biology: Regime Shifts in Ecology and Evolution"

[[4]]$defined_type
[1] "presentation"

[[4]]$status
[1] "Public"

[[4]]$version
[1] 2

[[4]]$published_date
[1] "21:04, Nov 14, 2012"

[[4]]$description
[1] "<p>Slides from my exit seminar satisfying the completion requirements for a PhD in Population Biology at UC Davis, November 13, 2012.</p>"

[[4]]$description_nohtml
[1] "Slides from my exit seminar satisfying the completion requirements for a PhD in Population Biology at UC Davis, November 13, 2012."

[[4]]$total_size
[1] "23.37 MB"

[[4]]$authors
[[4]]$authors[[1]]
[[4]]$authors[[1]]$first_name
[1] "Carl"

[[4]]$authors[[1]]$last_name
[1] "Boettiger"

[[4]]$authors[[1]]$id
[1] 96387

[[4]]$authors[[1]]$full_name
[1] "Carl Boettiger"



[[4]]$tags
[[4]]$tags[[1]]
[[4]]$tags[[1]]$id
[1] 49296

[[4]]$tags[[1]]$name
[1] " comparative methods"


[[4]]$tags[[2]]
[[4]]$tags[[2]]$id
[1] 49295

[[4]]$tags[[2]]$name
[1] "seminar"


[[4]]$tags[[3]]
[[4]]$tags[[3]]$id
[1] 49294

[[4]]$tags[[3]]$name
[1] " presentation"


[[4]]$tags[[4]]
[[4]]$tags[[4]]$id
[1] 49293

[[4]]$tags[[4]]$name
[1] " labrids"


[[4]]$tags[[5]]
[[4]]$tags[[5]]$id
[1] 49036

[[4]]$tags[[5]]$name
[1] " warning singals"


[[4]]$tags[[6]]
[[4]]$tags[[6]]$id
[1] 47771

[[4]]$tags[[6]]$name
[1] " phylogenetics"



[[4]]$categories
[[4]]$categories[[1]]
[[4]]$categories[[1]]$id
[1] 24

[[4]]$categories[[1]]$name
[1] "Evolutionary biology"


[[4]]$categories[[2]]
[[4]]$categories[[2]]$id
[1] 39

[[4]]$categories[[2]]$name
[1] "Ecology"



[[4]]$files
[[4]]$files[[1]]
[[4]]$files[[1]]$download_url
[1] "http://files.figshare.com/177712/boettiger.pdf"

[[4]]$files[[1]]$size
[1] "23.37 MB"

[[4]]$files[[1]]$id
[1] 177712

[[4]]$files[[1]]$mime_type
[1] "application/pdf"

[[4]]$files[[1]]$name
[1] "boettiger.pdf"



[[4]]$links
list()


[[5]]
[[5]]$article_id
[1] 97279

[[5]]$title
[1] "Presentation at the Computational Science Graduate Fellowship Conference (2012): Regime shifts in Ecology & Evolution"

[[5]]$defined_type
[1] "presentation"

[[5]]$status
[1] "Public"

[[5]]$version
[1] 2

[[5]]$published_date
[1] "00:19, Nov 09, 2012"

[[5]]$description
[1] "<p>Slides from my talk at CSGF 2012 in Washington DC. &nbsp;</p>\n<p>A video recording of the talk is available here:&nbsp;<a href=\"http://www.youtube.com/watch?v=xwIIVdyKe4o\">http://www.youtube.com/watch?v=xwIIVdyKe4o</a></p>"

[[5]]$description_nohtml
[1] "Slides from my talk at CSGF 2012 in Washington DC.A video recording of the talk is available here:http://www.youtube.com/watch?v=xwIIVdyKe4o"

[[5]]$total_size
[1] "13.39 MB"

[[5]]$authors
[[5]]$authors[[1]]
[[5]]$authors[[1]]$first_name
[1] "Carl"

[[5]]$authors[[1]]$last_name
[1] "Boettiger"

[[5]]$authors[[1]]$id
[1] 96387

[[5]]$authors[[1]]$full_name
[1] "Carl Boettiger"



[[5]]$tags
[[5]]$tags[[1]]
[[5]]$tags[[1]]$id
[1] 49084

[[5]]$tags[[1]]$name
[1] " hpc"


[[5]]$tags[[2]]
[[5]]$tags[[2]]$id
[1] 49035

[[5]]$tags[[2]]$name
[1] "regime shifts"


[[5]]$tags[[3]]
[[5]]$tags[[3]]$id
[1] 276

[[5]]$tags[[3]]$name
[1] "phylogenetics"


[[5]]$tags[[4]]
[[5]]$tags[[4]]$id
[1] 277

[[5]]$tags[[4]]$name
[1] "comparative methods"



[[5]]$categories
[[5]]$categories[[1]]
[[5]]$categories[[1]]$id
[1] 106

[[5]]$categories[[1]]$name
[1] "Science Policy"


[[5]]$categories[[2]]
[[5]]$categories[[2]]$id
[1] 24

[[5]]$categories[[2]]$name
[1] "Evolutionary biology"


[[5]]$categories[[3]]
[[5]]$categories[[3]]$id
[1] 34

[[5]]$categories[[3]]$name
[1] "Environmental science"


[[5]]$categories[[4]]
[[5]]$categories[[4]]$id
[1] 39

[[5]]$categories[[4]]$name
[1] "Ecology"


[[5]]$categories[[5]]
[[5]]$categories[[5]]$id
[1] 125

[[5]]$categories[[5]]$name
[1] "Bioinformatics"


[[5]]$categories[[6]]
[[5]]$categories[[6]]$id
[1] 77

[[5]]$categories[[6]]$name
[1] "Applied Computer Science"



[[5]]$files
[[5]]$files[[1]]
[[5]]$files[[1]]$download_url
[1] "http://files.figshare.com/143048/boettiger.pdf"

[[5]]$files[[1]]$size
[1] "13.39 MB"

[[5]]$files[[1]]$id
[1] 143048

[[5]]$files[[1]]$mime_type
[1] "application/pdf"

[[5]]$files[[1]]$name
[1] "boettiger.pdf"



[[5]]$links
[[5]]$links[[1]]
[[5]]$links[[1]]$link
[1] "http://www.youtube.com/watch?v=xwIIVdyKe4o"

[[5]]$links[[1]]$id
[1] 1255




[[6]]
[[6]]$article_id
[1] 97218

[[6]]$title
[1] "Regime shifts in ecology and evolution (PhD Dissertation)"

[[6]]$defined_type
[1] "paper"

[[6]]$status
[1] "Public"

[[6]]$version
[1] 1

[[6]]$published_date
[1] "18:34, Nov 06, 2012"

[[6]]$description
[1] "<p>The most pressing issues of our time are all characterized by sudden regime shifts: the collapse of&nbsp;marine fisheries or stock-markets, the overthrow of governments, shifts in global climate. Regime&nbsp;shifts, or sudden transitions in dynamical behavior of a system, underly many important phenomena in ecological and evolutionary problems. How do they arise? How can we identify when a&nbsp;shift has occurred? Can we forecast these shifts? Here I address each of these central questions in&nbsp;the context of a particular system. First, I show how stochasticity in eco-evolutionary dynamics&nbsp;can give rise two different domains, or regimes, governing the behavior of evolutionary trajectories (Boettiger et al., 2010). In the next chapter, I turn to the question of identifying evolutionary&nbsp;shifts from data using phylogenetic trees and morphological trait data of extant species (Boettiger&nbsp;et al., 2012). In the last chapter, I adapt the approach of the previous section which allowed me&nbsp;to quantify the information available in a given data set that could detect a shift into an approach&nbsp;for detecting regime shifts in ecological time series data before the occur (Boettiger and Hastings,&nbsp;2012).</p>\n<div>&nbsp;</div>"

[[6]]$description_nohtml
[1] "The most pressing issues of our time are all characterized by sudden regime shifts: the collapse of marine fisheries or stock-markets, the overthrow of governments, shifts in global climate. Regime shifts, or sudden transitions in dynamical behavior of a system, underly many important phenomena in ecological and evolutionary problems. How do they arise? How can we identify when a shift has occurred? Can we forecast these shifts? Here I address each of these central questions in the context of a particular system. First, I show how stochasticity in eco-evolutionary dynamics can give rise two different domains, or regimes, governing the behavior of evolutionary trajectories (Boettiger et al., 2010). In the next chapter, I turn to the question of identifying evolutionary shifts from data using phylogenetic trees and morphological trait data of extant species (Boettiger et al., 2012). In the last chapter, I adapt the approach of the previous section which allowed me to quantify the information available in a given data set that could detect a shift into an approach for detecting regime shifts in ecological time series data before the occur (Boettiger and Hastings, 2012)."

[[6]]$total_size
[1] "1.91 MB"

[[6]]$authors
[[6]]$authors[[1]]
[[6]]$authors[[1]]$first_name
[1] "Carl"

[[6]]$authors[[1]]$last_name
[1] "Boettiger"

[[6]]$authors[[1]]$id
[1] 96387

[[6]]$authors[[1]]$full_name
[1] "Carl Boettiger"



[[6]]$tags
[[6]]$tags[[1]]
[[6]]$tags[[1]]$id
[1] 49036

[[6]]$tags[[1]]$name
[1] " warning singals"


[[6]]$tags[[2]]
[[6]]$tags[[2]]$id
[1] 49035

[[6]]$tags[[2]]$name
[1] "regime shifts"



[[6]]$categories
[[6]]$categories[[1]]
[[6]]$categories[[1]]$id
[1] 24

[[6]]$categories[[1]]$name
[1] "Evolutionary biology"


[[6]]$categories[[2]]
[[6]]$categories[[2]]$id
[1] 39

[[6]]$categories[[2]]$name
[1] "Ecology"



[[6]]$files
[[6]]$files[[1]]
[[6]]$files[[1]]$download_url
[1] "http://files.figshare.com/142977/dissertation.pdf"

[[6]]$files[[1]]$size
[1] "1.91 MB"

[[6]]$files[[1]]$id
[1] 142977

[[6]]$files[[1]]$mime_type
[1] "application/pdf"

[[6]]$files[[1]]$name
[1] "dissertation.pdf"



[[6]]$links
[[6]]$links[[1]]
[[6]]$links[[1]]$link
[1] "http://carlboettiger.info"

[[6]]$links[[1]]$id
[1] 1235




[[7]]
[[7]]$article_id
[1] 96919

[[7]]$title
[1] "Lab Notebook, 2011"

[[7]]$defined_type
[1] "fileset"

[[7]]$status
[1] "Public"

[[7]]$version
[1] 1

[[7]]$published_date
[1] "22:44, Oct 29, 2012"

[[7]]$description
[1] "<p>Permanent archive of Carl Boettiger's open lab notebook entries for the year 2011 (<a href=\"http://www.carlboettiger.info/archives.html\">http://www.carlboettiger.info/archives.html</a>).&nbsp;Entries are archived in plain text UTF-8. Written in pandoc-flavored Markdown.&nbsp;Meets the goals of the Data Management Plan:&nbsp;<a href=\"http://www.carlboettiger.info/2012/10/09/data-management-plan.html\">http://www.carlboettiger.info/2012/10/09/data-management-plan.html</a></p>"

[[7]]$description_nohtml
[1] "Permanent archive of Carl Boettiger's open lab notebook entries for the year 2011 (http://www.carlboettiger.info/archives.htmlhttp://www.carlboettiger.info/2012/10/09/data-management-plan.html"

[[7]]$total_size
[1] "334.22 KB"

[[7]]$authors
[[7]]$authors[[1]]
[[7]]$authors[[1]]$first_name
[1] "Carl"

[[7]]$authors[[1]]$last_name
[1] "Boettiger"

[[7]]$authors[[1]]$id
[1] 96387

[[7]]$authors[[1]]$full_name
[1] "Carl Boettiger"



[[7]]$tags
[[7]]$tags[[1]]
[[7]]$tags[[1]]$id
[1] 48404

[[7]]$tags[[1]]$name
[1] "Ecology"


[[7]]$tags[[2]]
[[7]]$tags[[2]]$id
[1] 11917

[[7]]$tags[[2]]$name
[1] "ecology"


[[7]]$tags[[3]]
[[7]]$tags[[3]]$id
[1] 46723

[[7]]$tags[[3]]$name
[1] " open science"


[[7]]$tags[[4]]
[[7]]$tags[[4]]$id
[1] 47395

[[7]]$tags[[4]]$name
[1] " evolution"



[[7]]$categories
[[7]]$categories[[1]]
[[7]]$categories[[1]]$id
[1] 24

[[7]]$categories[[1]]$name
[1] "Evolutionary biology"


[[7]]$categories[[2]]
[[7]]$categories[[2]]$id
[1] 39

[[7]]$categories[[2]]$name
[1] "Ecology"



[[7]]$files
[[7]]$files[[1]]
[[7]]$files[[1]]$download_url
[1] "http://files.figshare.com/101976/2011.tar.gz"

[[7]]$files[[1]]$size
[1] "342 KB"

[[7]]$files[[1]]$id
[1] 101976

[[7]]$files[[1]]$mime_type
[1] "application/x-gzip"

[[7]]$files[[1]]$name
[1] "2011.tar.gz"



[[7]]$links
[[7]]$links[[1]]
[[7]]$links[[1]]$link
[1] "http://www.carlboettiger.info/2011"

[[7]]$links[[1]]$id
[1] 1148




[[8]]
[[8]]$article_id
[1] 96916

[[8]]$title
[1] "Lab Notebook, 2010"

[[8]]$defined_type
[1] "fileset"

[[8]]$status
[1] "Public"

[[8]]$version
[1] 1

[[8]]$published_date
[1] "22:39, Oct 29, 2012"

[[8]]$description
[1] "<p>Permanent archive of Carl Boettiger's open lab notebook entries for the year 2010 (<a href=\"http://www.carlboettiger.info/archives.html\">http://www.carlboettiger.info/archives.html</a>).&nbsp;Entries are archived in plain text UTF-8. Written in pandoc-flavored Markdown.&nbsp;Meets the goals of the Data Management Plan:&nbsp;<a href=\"http://www.carlboettiger.info/2012/10/09/data-management-plan.html\">http://www.carlboettiger.info/2012/10/09/data-management-plan.html</a></p>"

[[8]]$description_nohtml
[1] "Permanent archive of Carl Boettiger's open lab notebook entries for the year 2010 (http://www.carlboettiger.info/archives.htmlhttp://www.carlboettiger.info/2012/10/09/data-management-plan.html"

[[8]]$total_size
[1] "253.80 KB"

[[8]]$authors
[[8]]$authors[[1]]
[[8]]$authors[[1]]$first_name
[1] "Carl"

[[8]]$authors[[1]]$last_name
[1] "Boettiger"

[[8]]$authors[[1]]$id
[1] 96387

[[8]]$authors[[1]]$full_name
[1] "Carl Boettiger"



[[8]]$tags
[[8]]$tags[[1]]
[[8]]$tags[[1]]$id
[1] 48404

[[8]]$tags[[1]]$name
[1] "Ecology"


[[8]]$tags[[2]]
[[8]]$tags[[2]]$id
[1] 11917

[[8]]$tags[[2]]$name
[1] "ecology"


[[8]]$tags[[3]]
[[8]]$tags[[3]]$id
[1] 46723

[[8]]$tags[[3]]$name
[1] " open science"


[[8]]$tags[[4]]
[[8]]$tags[[4]]$id
[1] 47395

[[8]]$tags[[4]]$name
[1] " evolution"



[[8]]$categories
[[8]]$categories[[1]]
[[8]]$categories[[1]]$id
[1] 24

[[8]]$categories[[1]]$name
[1] "Evolutionary biology"


[[8]]$categories[[2]]
[[8]]$categories[[2]]$id
[1] 39

[[8]]$categories[[2]]$name
[1] "Ecology"



[[8]]$files
[[8]]$files[[1]]
[[8]]$files[[1]]$download_url
[1] "http://files.figshare.com/101975/2010.tar.gz"

[[8]]$files[[1]]$size
[1] "260 KB"

[[8]]$files[[1]]$id
[1] 101975

[[8]]$files[[1]]$mime_type
[1] "application/x-gzip"

[[8]]$files[[1]]$name
[1] "2010.tar.gz"



[[8]]$links
list()


[[9]]
[[9]]$article_id
[1] 95839

[[9]]$title
[1] "R code for the function fs_create "

[[9]]$defined_type
[1] "dataset"

[[9]]$status
[1] "Public"

[[9]]$version
[1] 1

[[9]]$published_date
[1] "19:48, Sep 13, 2012"

[[9]]$description
[1] "<p>An R implementation of the figshare API function to create a new article. &nbsp;Look Ma, I can share nice code on figshare!</p>\n<p>&nbsp;</p>\n<p>This function is part of the <a href=\"https://github.com/ropensci/rfigshare\">rfigshare</a> package by the <a href=\"http://ropensci.org\">rOpenSci</a> project. &nbsp;</p>"

[[9]]$description_nohtml
[1] "An R implementation of the figshare API function to create a new article.  Look Ma, I can share nice code on figshare!This function is part of therfigsharerOpenSci"

[[9]]$total_size
[1] "2.04 KB"

[[9]]$authors
[[9]]$authors[[1]]
[[9]]$authors[[1]]$first_name
[1] "Carl"

[[9]]$authors[[1]]$last_name
[1] "Boettiger"

[[9]]$authors[[1]]$id
[1] 96387

[[9]]$authors[[1]]$full_name
[1] "Carl Boettiger"



[[9]]$tags
[[9]]$tags[[1]]
[[9]]$tags[[1]]$id
[1] 47906

[[9]]$tags[[1]]$name
[1] " code"


[[9]]$tags[[2]]
[[9]]$tags[[2]]$id
[1] 47624

[[9]]$tags[[2]]$name
[1] "R"


[[9]]$tags[[3]]
[[9]]$tags[[3]]$id
[1] 42046

[[9]]$tags[[3]]$name
[1] "r"



[[9]]$categories
[[9]]$categories[[1]]
[[9]]$categories[[1]]$id
[1] 77

[[9]]$categories[[1]]$name
[1] "Applied Computer Science"



[[9]]$files
[[9]]$files[[1]]
[[9]]$files[[1]]$download_url
[1] "http://files.figshare.com/98377/fs_create.R"

[[9]]$files[[1]]$size
[1] "2 KB"

[[9]]$files[[1]]$id
[1] 98377

[[9]]$files[[1]]$mime_type
[1] "text/plain"

[[9]]$files[[1]]$name
[1] "fs_create.R"



[[9]]$links
list()


[[10]]
[[10]]$article_id
[1] 138

[[10]]$title
[1] "Labrid adaptive peaks"

[[10]]$defined_type
[1] "figure"

[[10]]$status
[1] "Public"

[[10]]$version
[1] 2

[[10]]$published_date
[1] "13:45, Dec 30, 2011"

[[10]]$description
[1] "<p>Described in the notebook: http://openwetware.org/wiki/User:Carl_Boettiger/Notebook/Comparative_Phylogenetics/2010/03/12</p>"

[[10]]$description_nohtml
[1] "Described in the notebook: http://openwetware.org/wiki/User:Carl_Boettiger/Notebook/Comparative_Phylogenetics/2010/03/12"

[[10]]$total_size
[1] "29.71 KB"

[[10]]$authors
[[10]]$authors[[1]]
[[10]]$authors[[1]]$first_name
[1] "Carl"

[[10]]$authors[[1]]$last_name
[1] "Boettiger"

[[10]]$authors[[1]]$id
[1] 96387

[[10]]$authors[[1]]$full_name
[1] "Carl Boettiger"



[[10]]$tags
[[10]]$tags[[1]]
[[10]]$tags[[1]]$id
[1] 277

[[10]]$tags[[1]]$name
[1] "comparative methods"


[[10]]$tags[[2]]
[[10]]$tags[[2]]$id
[1] 276

[[10]]$tags[[2]]$name
[1] "phylogenetics"


[[10]]$tags[[3]]
[[10]]$tags[[3]]$id
[1] 275

[[10]]$tags[[3]]$name
[1] "fins"


[[10]]$tags[[4]]
[[10]]$tags[[4]]$id
[1] 274

[[10]]$tags[[4]]$name
[1] "labrids"



[[10]]$categories
[[10]]$categories[[1]]
[[10]]$categories[[1]]$id
[1] 24

[[10]]$categories[[1]]$name
[1] "Evolutionary biology"


[[10]]$categories[[2]]
[[10]]$categories[[2]]$id
[1] 39

[[10]]$categories[[2]]$name
[1] "Ecology"



[[10]]$files
[[10]]$files[[1]]
[[10]]$files[[1]]$download_url
[1] "http://files.figshare.com/137/Labrid_fins.png"

[[10]]$files[[1]]$size
[1] "30 KB"

[[10]]$files[[1]]$id
[1] 137

[[10]]$files[[1]]$mime_type
[1] "image/png"

[[10]]$files[[1]]$name
[1] "Labrid_fins.png"



[[10]]$links
list()

```


Note that we can easily grab the ids with the wrapper function `fs_ids`:



```r
fs_ids(all_mine)
```

```
 [1] 105136 105135  97653  97500  97279  97218  96919  96916  95839    138
```




We can delete unwanted files that are not public with `fs_delete`:  


```r
fs_delete(id)
```



# Publishing On Figshare from R

Before you can use `rfigshare` effectively, you will need to set up your authentication credentials by obtaining a set of API keys from [FigShare.org](http://fishshare.org).  See our [Getting Started with rfigshare](https://github.com/ropensci/rfigshare/blob/master/inst/doc/getting_started.md) tutorial for a step by step guide.  


Now that you have created your credentials, we are ready to start posting content to FigShare using R.  

### Step 1: Create a new article


```r
require(rfigshare)
```



An article on FigShare can be a figure, poster, dataset, paper, or other media, and can contain an arbitrary number of files or attachements.  Each article has a unique ID number which we will use to interact with it.  All articles have the essential scientific metadata including a title, at least one author, and a description or abstract.  Categories can be selected from a fixed list, while articles can be assigned any tags.  This can be done step by step using dedicated functions, or simultaneously using the special function `fs_new_article`.  For example, the command


```r
id <- fs_new_article(title = "A Test of rfigshare", description = "This is a test of the fs_new_aricle function and related methods", 
    type = "figure", authors = c("Karthik Ram", "Scott Chamberlain"), tags = c("ecology", 
        "openscience"), categories = "Ecology", links = "http://ropensci.org", 
    files = "figure/rfigshare.png")
```

```
Your article has been created! Your id number is 95802
```

```
found ids for all authors
```


Creates a new article with the metadata given and returns the article id number, so we can make future modifications quickly.  

## Step 2: Examine and modify article

We can check out the details of our new article to confirm the successful creation:


```r
fs_details(id)
```

```
$article_id
[1] 95802

$title
[1] "A Test of rfigshare"

$views
[1] 0

$downloads
[1] 0

$shares
[1] 0

$doi
[1] "http://dx.doi.org/10.6084/m9.figshare.95802"

$defined_type
[1] "figure"

$status
[1] "Drafts"

$version
[1] 1

$published_date
[1] "23:18, Sep 10, 2012"

$description
[1] "This is a test of the fs_new_aricle function and related methods"

$total_size
[1] "17.78 KB"

$owner
$owner$id
[1] 96387

$owner$full_name
[1] "Carl Boettiger"


$authors
$authors[[1]]
$authors[[1]]$first_name
[1] "Carl"

$authors[[1]]$last_name
[1] "Boettiger"

$authors[[1]]$id
[1] 96387

$authors[[1]]$full_name
[1] "Carl Boettiger"


$authors[[2]]
$authors[[2]]$first_name
[1] "Karthik"

$authors[[2]]$last_name
[1] "Ram"

$authors[[2]]$id
[1] 97306

$authors[[2]]$full_name
[1] "Karthik Ram"


$authors[[3]]
$authors[[3]]$first_name
[1] "Scott"

$authors[[3]]$last_name
[1] "Chamberlain"

$authors[[3]]$id
[1] 96554

$authors[[3]]$full_name
[1] "Scott Chamberlain"



$tags
$tags[[1]]
$tags[[1]]$id
[1] 47864

$tags[[1]]$name
[1] "openscience"


$tags[[2]]
$tags[[2]]$id
[1] 11917

$tags[[2]]$name
[1] "ecology"



$categories
$categories[[1]]
$categories[[1]]$id
[1] 39

$categories[[1]]$name
[1] "Ecology"



$files
$files[[1]]
$files[[1]]$size
[1] "18 KB"

$files[[1]]$id
[1] 98324

$files[[1]]$mime_type
[1] "image/png"

$files[[1]]$name
[1] "rfigshare.png"



$links
$links[[1]]
$links[[1]]$link
[1] "http://ropensci.org"

$links[[1]]$id
[1] 673

```

Note that the submitter is automatically added as an author, though it will not hurt to specify them in the author list if you want anyway.  Note also the extra metadata we get, such as filesize and article views. (And hopefully we'll format that output to a pretty-printing version soon). 

If there is something we need to change, we can edit our article accordingly.  The `fs_update` function to modifies the title, description, and type, while `fs_add_tags`, `fs_add_categories`, `fs_add_authors` etc, can add missing data.  If we upload a new figure, it will overwrite this one.


```r
fs_update(id, title = "An awesome test of rfigshare")
```


# Step 3: Sharing or deleting your article

Once we are ready to share this, we can release the article privately or publicly.  We actually could have done this during our `fs_new_article` step by setting `visibility="private"`.  


```r
fs_make_private(id)
```

```
Response [http://api.figshare.com/v1/my_data/articles/95802/action/make_private]
  Status: 200
  Content-type: application/json; charset=UTF-8
{"success": "Article status changed to Private"} 
```


If we need to remove this example file we can delete it


```r
fs_delete(id)
```


Note that articles declared "public" cannot be deleted, and changes will appear as new versions of the same article.  `fs_details` can help you view particular versions of public articles.  
# Getting Started with rfigshare

![](http://farm9.staticflickr.com/8180/7950489358_ea902bdaae_o.png)


## Obtaining your API keys

Note that there is a nice video introduction to creating applications for the API on the [figshare blog](http://figshare.com/blog/figshare_API_available_to_all/48).  The following tutorial provides a simple walkthrough of how to go about getting your figshare API keys set up so that you can use the `rfigshare` package.  


Create a user account on [FigShare](http://figshare.com) and log in.  From your homepage, select "Applications" from the drop-down menu,

![](http://farm9.staticflickr.com/8171/7950489558_5172515057_o.png)

Create a new application:

![](http://farm9.staticflickr.com/8038/7950490158_7feaf680bd_o.png)


Enter in the following information: 

![](http://farm9.staticflickr.com/8305/7950490562_02846cea92_o.png)

Then navigate over to the permissions tab.  To get the most out of `rfigshare` you'll want to enable all permissions:

![](http://farm9.staticflickr.com/8448/7950491064_c3820e62bd_o.png)

Save the new settings, and then open the application again (View/Edit menu) and click on the "Access Codes" tab.

![](http://farm9.staticflickr.com/8308/7950491470_621da9c5d1_o.png)

Record each if the keys into R as follows.  You might want to put this bit of R code into your `.Rprofile` to avoid entering it each time in the future:

```r
options(FigshareKey = "qMDabXXXXXXXXXXXXXXXXX")
options(FigsharePrivateKey = "zQXXXXXXXXXXXXXXXXXXXX")
options(FigshareToken = "SrpxabQXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX")
options(FigsharePrivateToken = "yqXXXXXXXXXXXXXXXXXXXX")
```

That's it! You are now ready to start using figshare.  Recall you can install the package directly from Github using: 

```r
require(devtools)
install_github("rfigshare", "ropensci")
```

Try authenticating with your credentials:


```r
require(rfigshare)
```

```
## Loading required package: rfigshare
```

```r
fs_auth()
```



Try a search for an author, or get the details on a paper:


```r
fs_author_search("Boettiger")
```

```
## Response [http://api.figshare.com/v1/my_data/authors?search_for=Boettiger]
##   Status: 200
##   Content-type: application/json; charset=UTF-8
## {"pages": 0, "results": 1, "start": 0, "per_page": 10, "items": [{"id": "96387", "fname": "Carl", "lname": "Boettiger", "full_name": "Carl Boettiger", "job_title": "", "description": "", "facebook": "", "twitter": "", "active": 1}]} 
```

```r
fs_details("138")
```

```
## Warning: text_content() deprecated. Use parsed_content(x, as = 'parsed')
```

```
## Loading required package: rjson
```

```
## $article_id
## [1] 138
## 
## $title
## [1] "Labrid adaptive peaks"
## 
## $views
## [1] 56
## 
## $downloads
## [1] 0
## 
## $shares
## [1] 0
## 
## $doi
## [1] "http://dx.doi.org/10.6084/m9.figshare.138"
## 
## $defined_type
## [1] "figure"
## 
## $status
## [1] "Public"
## 
## $version
## [1] 1
## 
## $published_date
## [1] "13:45, Dec 30, 2011"
## 
## $description
## [1] "Described in the notebook: http://openwetware.org/wiki/User:Carl_Boettiger/Notebook/Comparative_Phylogenetics/2010/03/12"
## 
## $total_size
## [1] "29.71 KB"
## 
## $owner
## $owner$id
## [1] 96387
## 
## $owner$full_name
## [1] "Carl Boettiger"
## 
## 
## $authors
## $authors[[1]]
## $authors[[1]]$first_name
## [1] "Carl"
## 
## $authors[[1]]$last_name
## [1] "Boettiger"
## 
## $authors[[1]]$id
## [1] 96387
## 
## $authors[[1]]$full_name
## [1] "Carl Boettiger"
## 
## 
## 
## $tags
## $tags[[1]]
## $tags[[1]]$id
## [1] 277
## 
## $tags[[1]]$name
## [1] "comparative methods"
## 
## 
## $tags[[2]]
## $tags[[2]]$id
## [1] 276
## 
## $tags[[2]]$name
## [1] "phylogenetics"
## 
## 
## $tags[[3]]
## $tags[[3]]$id
## [1] 275
## 
## $tags[[3]]$name
## [1] "fins"
## 
## 
## $tags[[4]]
## $tags[[4]]$id
## [1] 274
## 
## $tags[[4]]$name
## [1] "labrids"
## 
## 
## 
## $categories
## $categories[[1]]
## $categories[[1]]$id
## [1] 24
## 
## $categories[[1]]$name
## [1] "Evolutionary biology"
## 
## 
## $categories[[2]]
## $categories[[2]]$id
## [1] 39
## 
## $categories[[2]]$name
## [1] "Ecology"
## 
## 
## 
## $files
## $files[[1]]
## $files[[1]]$size
## [1] "30 KB"
## 
## $files[[1]]$id
## [1] 137
## 
## $files[[1]]$mime_type
## [1] "image/png"
## 
## $files[[1]]$name
## [1] "Labrid_fins.png"
## 
## 
## 
## $links
## list()
## 
```


Try creating your own content:


```r
fs_create("Test title", "description of test", "dataset")
```

```
## Warning: text_content() deprecated. Use parsed_content(x, as = 'parsed')
```

```
## Your article has been created! Your id number is 95717
```

```
## [1] 95717
```


This creates an article with the essential metadata information.  In the next tutorial, [Publishing on FigShare from R](https://github.com/ropensci/rfigshare/blob/master/inst/doc/publishing_on_figshare.md) we will describe how to add files, tags, categories, authors, and links to your draft, and then publish it either privately or publicly.   





---
title: rfigshare tutorial
author: Carl Boettiger
---


<!--
%\VignetteEngine{knitr::knitr}
%\VignetteIndexEntry{An Introduction to the rfigshare package}
-->

[![Build Status](https://api.travis-ci.org/ropensci/rfigshare.png)](https://travis-ci.org/ropensci/rfigshare)




rfigshare
==========

*An R interface to [FigShare](http://figshare.com)*

* Maintainer: Carl Boettiger, [cboettig](https://github.com/cboettig)
* License: [CC0](http://creativecommons.org/publicdomain/zero/1.0/)
* Contact: Report bugs, questions, or feature requests on the [Issues Tracker](https://github.com/ropensci/rfigshare/issues), or get in touch with us at [info@ropensci.org](mailto:info@ropensci.org)


Getting Started
---------------
Figshare is an online digital repository where researchers can preserve and share their research outputs, including figures, datasets, images, and videos. It is free to upload content and free to access, in adherence to the principle of open data.

Key Features:
- Showcase your institution's research with a customizable portal of all public research outputs using the reporting and statistics feature.
- Have full control of your institution's research outputs with private storage, public storage and collaborative spaces with the data management feature.
- Filter your institution's research by department, category or file type, and rank content by most viewed, downloaded or shared with the data dissemination feauture. 
- Manage the curation of files to be made public, control quotas, and administer rights with the user group administration feature. 





# Using rfigshare


```{r}
library("rfigshare")
```


```{r include = FALSE}
# This loads the rOpenSci figshare sandbox credentials, so that the example 
# can run automatically during check and install.  Unlike normal figshare accounts,
# data loaded to this testing sandbox is periodically purged.  
fs_auth(token = "xdBjcKOiunwjiovwkfTF2QjGhROeLMw0y0nSCSgvg3YQxdBjcKOiunwjiovwkfTF2Q", token_secret = "4mdM3pfekNGO16X4hsvZdg")
```

The first time you use an `rfigshare` function, it will ask you to authenticate online. Just log in and click okay to authenticate rfigshare.  R will allow you to cache your login credentials so that you won't be asked to authenticate again (even between R sessions), as long as you are using the same working directory in future.  

Try a search for an author:


```{r}
fs_author_search("Boettiger")
```



Try creating your own content:


```{r}
id <- fs_create("Test title", "description of test")
```


This creates an article with the essential metadata information. We can now upload the dataset to this newly created figshare object using `fs_upload`.  For the purposes of this example, we'll just upload one of R's built-in datasets:


```{r}
data(mtcars)
write.csv(mtcars, "mtcars.csv")
fs_upload(id, "mtcars.csv")
```


Note that we must pass the id number of our the newly created figshare object as the first argument.  Similar functions let us add additional authors, tags, categories, and links, e.g.


```{r}
fs_add_tags(id, "demo")
```



Minimal metadata includes title, description, type, and at least one tag and one category.  We can add categories using either the category id or the name of the category, but it must be one of the preset categories available.  We can ask the API for a list of all the categories:


```{r}
fs_category_list()
```


And we can add the category or categories we like,


```{r}
fs_add_categories(id, c("Education", "Software Engineering"))
```



The file we have created remains saved as a draft until we publish it, either publicly or privately.  Note that once a file is declared public, it cannot be made private or deleted.  Let's release this dataset privately:


```{r}
fs_make_private(id)
```


We can now share the dataset with collaborators by way of the private url.  

### The quick and easy way

The `rfigshare` package will let you create a new figshare article with additional authors, tags, categories, etc in a single command usnig the `fs_new_article` function.  The essential metadata `title`, `description` and `type` are required, but any other information we omit at this stage can be added later.  If we set `visibility` to private or public, the article is published on figshare immediately.  


```{r}
data(mtcars)
write.csv(mtcars,"mtcars.csv")
id <- fs_new_article(title="A Test of rfigshare", 
                     description="This is a test", 
                     type="dataset", 
                     authors=c("Karthik Ram", "Scott Chamberlain"), 
                     tags=c("ecology", "openscience"), 
                     categories="Ecology", 
                     links="http://ropensci.org", 
                     files="mtcars.csv",
                     visibility="private")
unlink("mtcars.csv") # clean up
```


# Examining Data on Figshare

We can view all available metadata of a figshare object. 


```{r}
fs_details(id)
```

You can see all of the files you have (Currently only up to 10):


```{r}
mine <- fs_browse()
mine[1:2]
```

Note that we can easily grab the ids with the wrapper function `fs_ids`:



```{r}
fs_ids(mine)
```


We can delete unwanted files that are not public with `fs_delete`:  


```{r}
fs_delete(id)
```

To cite package `rfigshare` in publications use:


```{r}
citation("rfigshare")
```


[![ropensci.org logo](http://ropensci.org/public_images/github_footer.png)](http://ropensci.org)\
<!--
%\VignetteEngine{knitr}
%\VignetteIndexEntry{A Markdown Vignette with knitr}
-->


# rfigshare - Share and explore figures, data, and publications on FigShare using R


![](http://farm9.staticflickr.com/8180/7950489358_ea902bdaae_o.png)


```{r include=FALSE}
opts_chunk$set(comment=NA, tidy=FALSE, warning=FALSE)
```



# Obtaining your API keys

There is a nice video introduction to creating applications for the API on the [figshare blog](http://figshare.com/blog/figshare_API_available_to_all/48).  The following tutorial provides a simple walkthrough of how to go about getting your figshare API keys set up so that you can use the `rfigshare` package.  


Create a user account on [FigShare](http://figshare.com) and log in.  From your homepage, select "Applications" from the drop-down menu,

![](http://farm9.staticflickr.com/8171/7950489558_5172515057_o.png)

Create a new application:

![](http://farm9.staticflickr.com/8038/7950490158_7feaf680bd_o.png)


Enter in the following information: 

![](http://farm9.staticflickr.com/8305/7950490562_02846cea92_o.png)

Then navigate over to the permissions tab.  To get the most out of `rfigshare` you'll want to enable all permissions:

![](http://farm9.staticflickr.com/8448/7950491064_c3820e62bd_o.png)

Save the new settings, and then open the application again (View/Edit menu) and click on the "Access Codes" tab.

![](http://farm9.staticflickr.com/8308/7950491470_621da9c5d1_o.png)

Record each if the keys into R as follows.  You might want to put this bit of R code into your `.Rprofile` to avoid entering it each time in the future:

```r
options(FigshareKey = "qMDabXXXXXXXXXXXXXXXXX")
options(FigsharePrivateKey = "zQXXXXXXXXXXXXXXXXXXXX")
options(FigshareToken = "SrpxabQXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX")
options(FigsharePrivateToken = "yqXXXXXXXXXXXXXXXXXXXX")
```


## Installing the R package

Now that we have the API credentials in place, actually using `rfigshare` is pretty easy.  Install the latest version of package directly from Github using: 

```r
require(devtools)
install_github("rfigshare", "ropensci")
```

# Using rfigshare


Try authenticating with your credentials:

``` {r }
require(rfigshare)
fs_auth()
````


Try a search for an author:

``` {r }
fs_author_search("Boettiger")
````

Try creating your own content:

``` {r }
id <- fs_create("Test title", "description of test", "dataset")
````

This creates an article with the essential metadata information. We can now upload the dataset to this newly created figshare object using `fs_upload`.  For the purposes of this example, we'll just upload one of R's built-in datasets:

``` {r }
data(mtcars)
write.csv(mtcars, "mtcars.csv")

fs_upload(id, "mtcars.csv")
```

Not that we must pass the id number of our the newly created figshare object as the first argument.  Similar functions let us add additional authors, tags, categories, and links, e.g.

``` {r }
fs_add_tags(id, "demo")
```


Minimal metadata includes title, description, type, and at least one tag and one category.  We can add categories using either the category id or the name of the category, but it must be one of the preset categories available.  We can ask the API for a list of all the categories:

``` {r results="hide"}
fs_category_list()
```

``` {r results="asis", echo=FALSE}
print(xtable::xtable(fs_category_list()), type="html")
```

And we can add the category or categories we like,

``` {r }
fs_add_categories(id, c("Education", "Software Engineering"))
```


The file we have created remains saved as a draft until we publish it, either publicly or privately.  Note that once a file is declared public, it cannot be made private or deleted.  Let's release this dataset privately:

``` {r }
fs_make_private(id)
```

We can now share the dataset with collaborators by way of the private url.  

### The quick and easy way

The `rfigshare` package will let you create a new figshare article with additional authors, tags, categories, etc in a single command usnig the `fs_new_article` function.  The essential metadata `title`, `description` and `type` are required, but any other information we omit at this stage can be added later.  If we set `visibility` to private or public, the article is published on figshare immediately.  

``` {r tidy=FALSE}
id <- fs_new_article(title="A Test of rfigshare", 
                     description="This is a test", 
                     type="figure", 
                     authors=c("Karthik Ram", "Scott Chamberlain"), 
                     tags=c("ecology", "openscience"), 
                     categories="Ecology", 
                     links="http://ropensci.org", 
                     files="figure/rfigshare.png",
                     visibility="private")
```

# Examining Data on Figshare

We can view all available metadata of a figshare object.  If the object is not public (e.g. draft or private), we have to add the `mine=TRUE` option

``` {r }
fs_details(id, mine=TRUE)
```

You can see all of the files you have:

``` {r } 
fs_browse(mine=TRUE)
```

Note that we can easily grab the ids with the wrapper function `fs_ids`:


``` {r }
fs_ids(all_mine)
```



We can delete unwanted files that are not public with `fs_delete`:  

``` {r }
fs_delete(id)
```
<!--
%\VignetteEngine{knitr}
%\VignetteIndexEntry{A Markdown Vignette with knitr}
-->


rOpenSci, Figshare and knitr.
========================================================

This is a tutorial of how you can create reproducible documents using [knitr]( http://yihui.name/knitr/) and pandoc, and seamlessly upload them to [figshare](http;//www.figshare.com) attaching a citable DOI to them. This document will walk you through the process of creating a document in knitr and uploading a compiled PDF to Figshare.  I make the following assumptions about your knowledge:

* You have set-up an account at [Figshare.com](http://www.figshare.com)

* You have installed [rfigshare](http://cran.r-project.org/web/packages/rfigshare/index.html), the [knitr](http://yihui.name/knitr/) package and are familiar with the concepts of knitr or sweave, as well as [pandoc](http://johnmacfarlane.net/pandoc/), for conversion to pdf (although this is only necessary if you want to convert your document to a pdf)

* You have successfully set-up the your credentials for rfigshare.  If not go to our [tutorial](http://github.com/ropensci/rfigshare/blob/master/inst/doc/tutorial.md) and make sure your credentials are properly set.


The goal of this document is to demonstrate how one could carry out a project using tools from [rOpenSci](http://www.ropensci.org/), knitr, and share the results on figshare in one continuous workflow.  To do this I'll be using a tutorial from one of our packages, [treebase](http://cran.r-project.org/web/packages/treebase/), which allows you to download trees from [TreebaseWEB](http://treebase.org/treebase-web/home.html;jsessionid=A258F89FBF584F44E0CDB740B8ECF3A8)


First I'll turn the cache on.
```{r headerchunk }
opts_chunk$set(cache=TRUE, autodep=TRUE)
dep_auto()

```
Then you'll want to download some data, and maybe make a plot, and say some things about how great your plot is.
```{r messages=FALSE,warning=FALSE,message=FALSE,fig.cap="My amazing tree!",fig.width=7,fig.height=4}
library(treebase)
tree <- search_treebase("Derryberry", "author")[[1]] 
#plotting only part of the tree because it's so large
plot.phylo(tree,y.lim=c(0,20))

```

Once you've made all your plots, and said all you want to say it's time to convert your document, and then create a new article using `fs_new_article()`


```{r eval = FALSE}
library(knitr)
library(rfigshare)
options(FigshareKey = "XXXXXXXX")
options(FigsharePrivateKey = "XXXXXXXX")
options(FigshareToken = "XXXXXXXX")
options(FigsharePrivateToken = "XXXXXXXX")

#knit document to pandoc markdown
knit("rfigtutorial.Rmd")
#convert to pdf
system("pandoc -S rfigtutorial.md -o rfigtutorial.pdf")

id <- fs_new_article(title="An rfigshare tutorial", 
                      description="How to create a document in knitr and 
                      upload it to figshare.com", 
                     
                      type="paper", 
                      authors=c("Edmund Hart"), 
                      tags=c("ecology", "openscience"), 
                      categories="Ecology", 
                      links="http://emhart.github.com", 
                      files="rfigtutorial.pdf",
                      visibility="draft")
```

The main advantage of this approach is that manuscripts can be worked on from within the R environment and then seemlessly uploaded to figshare. Also it's best practice to store your key values in your `.Rprofile` so I would reccomend file Be sure to run `fs_make_public(id)` when you're ready to make your article public.  


<!--
%\VignetteEngine{knitr}
%\VignetteIndexEntry{A Markdown Vignette with knitr}
-->


`ro warning=FALSE, comment=NA or`

# Publishing On Figshare from R

Before you can use `rfigshare` effectively, you will need to set up your authentication credentials by obtaining a set of API keys from [FigShare.org](http://fishshare.org).  See our [Getting Started with rfigshare](https://github.com/ropensci/rfigshare/blob/master/inst/doc/getting_started.md) tutorial for a step by step guide.  


Now that you have created your credentials, we are ready to start posting content to FigShare using R.  

### Step 1: Create a new article

``` {r }
require(rfigshare)
```


An article on FigShare can be a figure, poster, dataset, paper, or other media, and can contain an arbitrary number of files or attachements.  Each article has a unique ID number which we will use to interact with it.  All articles have the essential scientific metadata including a title, at least one author, and a description or abstract.  Categories can be selected from a fixed list, while articles can be assigned any tags.  This can be done step by step using dedicated functions, or simultaneously using the special function `fs_new_article`.  For example, the command

``` {r }
id <- fs_new_article(title="A Test of rfigshare", 
                     description="This is a test of the fs_new_aricle function and related methods", 
                     type="figure", 
                     authors=c("Karthik Ram", "Scott Chamberlain"), 
                     tags=c("ecology", "openscience"), 
                     categories="Ecology", 
                     links="http://ropensci.org", 
                     files="figure/rfigshare.png")
```

Creates a new article with the metadata given and returns the article id number, so we can make future modifications quickly.  

## Step 2: Examine and modify article

We can check out the details of our new article to confirm the successful creation:

``` {r }
fs_details(id)

```
Note that the submitter is automatically added as an author, though it will not hurt to specify them in the author list if you want anyway.  Note also the extra metadata we get, such as filesize and article views. (And hopefully we'll format that output to a pretty-printing version soon). 

If there is something we need to change, we can edit our article accordingly.  The `fs_update` function to modifies the title, description, and type, while `fs_add_tags`, `fs_add_categories`, `fs_add_authors` etc, can add missing data.  If we upload a new figure, it will overwrite this one.

```{r }
fs_update(id, title="An awesome test of rfigshare")
```

# Step 3: Sharing or deleting your article

Once we are ready to share this, we can release the article privately or publicly.  We actually could have done this during our `fs_new_article` step by setting `visibility="private"`.  

``` {r }
fs_make_private(id)
```

If we need to remove this example file we can delete it

``` {r }
fs_delete(id)
```

Note that articles declared "public" cannot be deleted, and changes will appear as new versions of the same article.  `fs_details` can help you view particular versions of public articles.  
<!--
%\VignetteEngine{knitr}
%\VignetteIndexEntry{A Markdown Vignette with knitr}
-->


# Getting Started with rfigshare

![](http://farm9.staticflickr.com/8180/7950489358_ea902bdaae_o.png)


## Obtaining your API keys

Note that there is a nice video introduction to creating applications for the API on the [figshare blog](http://figshare.com/blog/figshare_API_available_to_all/48).  The following tutorial provides a simple walkthrough of how to go about getting your figshare API keys set up so that you can use the `rfigshare` package.  


Create a user account on [FigShare](http://figshare.com) and log in.  From your homepage, select "Applications" from the drop-down menu,

![](http://farm9.staticflickr.com/8171/7950489558_5172515057_o.png)

Create a new application:

![](http://farm9.staticflickr.com/8038/7950490158_7feaf680bd_o.png)


Enter in the following information: 

![](http://farm9.staticflickr.com/8305/7950490562_02846cea92_o.png)

Then navigate over to the permissions tab.  To get the most out of `rfigshare` you'll want to enable all permissions:

![](http://farm9.staticflickr.com/8448/7950491064_c3820e62bd_o.png)

Save the new settings, and then open the application again (View/Edit menu) and click on the "Access Codes" tab.

![](http://farm9.staticflickr.com/8308/7950491470_621da9c5d1_o.png)

Record each if the keys into R as follows.  You might want to put this bit of R code into your `.Rprofile` to avoid entering it each time in the future:

```r
options(FigshareKey = "qMDabXXXXXXXXXXXXXXXXX")
options(FigsharePrivateKey = "zQXXXXXXXXXXXXXXXXXXXX")
options(FigshareToken = "SrpxabQXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX")
options(FigsharePrivateToken = "yqXXXXXXXXXXXXXXXXXXXX")
```

That's it! You are now ready to start using figshare.  Recall you can install the package directly from Github using: 

```r
require(devtools)
install_github("rfigshare", "ropensci")
```

Try authenticating with your credentials:

``` {r }
require(rfigshare)
fs_auth()
````


Try a search for an author, or get the details on a paper:

``` {r }
fs_author_search("Boettiger")
fs_details("138")
````

Try creating your own content:

``` {r }
fs_create("Test title", "description of test", "dataset")
````

This creates an article with the essential metadata information.  In the next tutorial, [Publishing on FigShare from R](https://github.com/ropensci/rfigshare/blob/master/inst/doc/publishing_on_figshare.md) we will describe how to add files, tags, categories, authors, and links to your draft, and then publish it either privately or publicly.   





<!--
%\VignetteEngine{knitr::knitr}
%\VignetteIndexEntry{An Introduction to the rfigshare package}
-->

[![Build Status](https://api.travis-ci.org/ropensci/rfigshare.png)](https://travis-ci.org/ropensci/rfigshare)




rfigshare
==========

*An R interface to [FigShare](http://figshare.com)*

* Maintainer: Carl Boettiger, [cboettig](https://github.com/cboettig)
* License: [CC0](http://creativecommons.org/publicdomain/zero/1.0/)
* Contact: Report bugs, questions, or feature requests on the [Issues Tracker](https://github.com/ropensci/rfigshare/issues), or get in touch with us at [info@ropensci.org](mailto:info@ropensci.org)

Installation
------------


```{r eval=FALSE}
require(devtools)
install_github("rfigshare", "ropensci")
```

Getting Started
---------------



# Using rfigshare


```{r}
require(rfigshare)
```


```{r include = FALSE}
# This loads the rOpenSci figshare sandbox credentials, so that the example 
# can run automatically during check and install.  Unlike normal figshare accounts,
# data loaded to this testing sandbox is periodically purged.  
fs_auth(token = "xdBjcKOiunwjiovwkfTF2QjGhROeLMw0y0nSCSgvg3YQxdBjcKOiunwjiovwkfTF2Q", token_secret = "4mdM3pfekNGO16X4hsvZdg")
```

The first time you use an `rfigshare` function, it will ask you to authenticate online. Just log in and click okay to authenticate rfigshare.  R will allow you to cache your login credentials so that you won't be asked to authenticate again (even between R sessions), as long as you are using the same working directory in future.  

Try a search for an author:


```{r}
fs_author_search("Boettiger")
```



Try creating your own content:


```{r}
id <- fs_create("Test title", "description of test")
```


This creates an article with the essential metadata information. We can now upload the dataset to this newly created figshare object using `fs_upload`.  For the purposes of this example, we'll just upload one of R's built-in datasets:


```{r}
data(mtcars)
write.csv(mtcars, "mtcars.csv")
fs_upload(id, "mtcars.csv")
```


Not that we must pass the id number of our the newly created figshare object as the first argument.  Similar functions let us add additional authors, tags, categories, and links, e.g.


```{r}
fs_add_tags(id, "demo")
```



Minimal metadata includes title, description, type, and at least one tag and one category.  We can add categories using either the category id or the name of the category, but it must be one of the preset categories available.  We can ask the API for a list of all the categories:


```{r}
fs_category_list()
```


And we can add the category or categories we like,


```{r}
fs_add_categories(id, c("Education", "Software Engineering"))
```



The file we have created remains saved as a draft until we publish it, either publicly or privately.  Note that once a file is declared public, it cannot be made private or deleted.  Let's release this dataset privately:


```{r}
fs_make_private(id)
```


We can now share the dataset with collaborators by way of the private url.  

### The quick and easy way

The `rfigshare` package will let you create a new figshare article with additional authors, tags, categories, etc in a single command usnig the `fs_new_article` function.  The essential metadata `title`, `description` and `type` are required, but any other information we omit at this stage can be added later.  If we set `visibility` to private or public, the article is published on figshare immediately.  


```{r}
data(mtcars)
write.csv(mtcars,"mtcars.csv")
id <- fs_new_article(title="A Test of rfigshare", 
                     description="This is a test", 
                     type="dataset", 
                     authors=c("Karthik Ram", "Scott Chamberlain"), 
                     tags=c("ecology", "openscience"), 
                     categories="Ecology", 
                     links="http://ropensci.org", 
                     files="mtcars.csv",
                     visibility="private")
unlink("mtcars.csv") # clean up
```


# Examining Data on Figshare

We can view all available metadata of a figshare object. 


```{r}
fs_details(id)
```

You can see all of the files you have (Currently only up to 10):


```{r}
mine <- fs_browse()
mine[1:2]
```

Note that we can easily grab the ids with the wrapper function `fs_ids`:



```{r}
fs_ids(mine)
```

```
 [1] 105136 105135  97653  97500  97279  97218  96919  96916  95839    138
```




We can delete unwanted files that are not public with `fs_delete`:  


```{r}
fs_delete(id)
```

To cite package `rfigshare` in publications use:


```{r}
citation("rfigshare")
```


[![](http://ropensci.org/public_images/github_footer.png)](http://ropensci.org)
---
title: "Updating public repositories"
author: "Martin John Hadley, @martinjhnhadley"
date: "`r Sys.Date()`"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Vignette Title}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

The process for updating an existing public repository with new versions of existing files requires multiple steps:

1. Delete existing copies of files
2. Upload new copies of the files
3. Make your changes public

Note that the example below is for a specific repository, as you are not the author of the article you will not be able to run the code without errors.

## Delete existing copies of files

Obtain the article_id of the deposit, this is the numeric component of the DOI after figshare e.g. https://doi.org/10.6084/m9.figshare.3761562

```{r, eval=FALSE, echo=TRUE}
library(rfigshare)
article_id <- 3761562
deposit_details <- fs_details(article_id)
deposit_details$title
```

Several files in this deposit are updated nightly, for instance:

```{r, eval=FALSE, echo=TRUE}
"OLIdata_YYYY-MM-DD.txt"
```

To delete this file, the file_id must be found. It is simplest to convert the lists to a `data_frame` such that they may be operated on with `dplyr`.

```{r, eval=FALSE, echo=TRUE}
library(dplyr)
deposit_files <- unlist(deposit_details$files)
deposit_files <- data.frame(split(deposit_files, names(deposit_files)),stringsAsFactors = F)
file_id <- deposit_files %>%
  filter(grepl("^OLIdata_", name)) %>%
  select(id) %>%
  .[[1]]
```

Prepare the article for the new version of the file, by deleting the existant version with `fs_delete`

```{r, eval=FALSE, echo=TRUE}
fs_delete(article_id, file_id)
```

## Upload new version of the file

The new file can be downloaded as follows:

```{r, eval=FALSE, echo=TRUE}
# This file does not exist in these training materials.
fs_upload(article_id, paste0("OLIdata_", as.Date(Sys.time())))
```

## Make changes public

The actions you have performed have been saved as draft changes, you must use `fs_make_public` to update the article and create a new version:

```{r, eval=FALSE, echo=TRUE}
fs_make_public(article_id)
```




% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fs_upload.R
\name{fs_upload_one}
\alias{fs_upload_one}
\title{Upload file to an article}
\usage{
fs_upload_one(article_id, file, session = fs_get_auth())
}
\arguments{
\item{article_id}{number}

\item{file}{path to file to upload}

\item{session}{the authentication credentials from \code{\link{fs_auth}} (optional)}
}
\description{
Upload file to an article
}
\details{
Article must be a draft, i.e. created by \code{\link{fs_create}} and not yet made public or private. Only articles of type "fileset" can have multiple files uploaded.
}
\examples{
\dontrun{
id <- fs_create("Title", "description", "figure")
fs_upload(id, "file.png")
} 
}
\references{
\url{http://api.figshare.com}
}
\seealso{
\code{\link{fs_auth}}
}
\author{
Carl Boettiger \email{cboettig@gmail.com}
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fs_create_author.R
\name{fs_create_author}
\alias{fs_create_author}
\title{Creates a figshare author}
\usage{
fs_create_author(full_name, session = fs_get_auth(), debug = FALSE)
}
\arguments{
\item{full_name}{full name of the author to create}

\item{session}{(optional) the authentication credentials from \code{\link{fs_auth}}. If not provided, will attempt to load from cache as long as figshare_auth has been run.}

\item{debug}{return PUT request visibly?}
}
\value{
author ID numbers
}
\description{
Creates a figshare author
}
\examples{
\dontrun{
fs_create_author("Benjamin Franklin") 
} 
}
\references{
\url{http://api.figshare.com}
}
\seealso{
\code{\link{fs_auth}}
}
\author{
Carl Boettiger \email{cboettig@gmail.com}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fs_make_public.R
\name{fs_make_public}
\alias{fs_make_public}
\title{Make an article public (for private or draft articles)}
\usage{
fs_make_public(article_id, session = fs_get_auth())
}
\arguments{
\item{article_id}{the id number of the article}

\item{session}{(optional) the authentication credentials from \code{\link{fs_auth}}.}
}
\value{
output of PUT request (invisibly)
}
\description{
Make an article public (for private or draft articles)
}
\details{
This function will make a draft or private article public, assigning it a DOI and making it permanently available through Figshare. If you use \code{\link{fs_upload}} to add new files to an existing public deposit, you must then use \code{fs_make_public} for those changes to be made in the public version of the repository.
}
\note{
NOTE: Public articles are assigned DOIs and cannot be deleted or made private once declared public! Public articles do not count against your quota space.
}
\examples{
\dontrun{
fs_make_public(123)
}
}
\references{
\url{http://api.figshare.com}
}
\seealso{
\code{\link{fs_auth}}
}
\author{
Carl Boettiger \email{cboettig@gmail.com}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fs_ids.R
\name{fs_ids}
\alias{fs_ids}
\title{Get a list of article id numbers from a search return}
\usage{
fs_ids(object)
}
\arguments{
\item{object}{the output of a search}
}
\value{
a list of article id numbers
}
\description{
Get a list of article id numbers from a search return
}
\examples{
\dontrun{
figshare_category() 
}
}
\references{
\url{http://api.figshare.com}
}
\author{
Carl Boettiger \email{cboettig@gmail.com}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fs_search.R
\name{fs_search}
\alias{fs_search}
\title{Advanced Search.}
\usage{
fs_search(query, author = NA, title = NA, description = NA,
  tag = NA, category = NA, from_date = NA, to_date = NA,
  mine = FALSE, public_only = FALSE, private_only = FALSE,
  drafts_only = FALSE, session = fs_get_auth(),
  base = "http://api.figshare.com/v1", debug = FALSE)
}
\arguments{
\item{query}{the search query}

\item{author}{Show only results by this author}

\item{title}{Show only results matching or partially matching this title}

\item{description}{Show only results matching or partially matching this description}

\item{tag}{Show only results matching this tag}

\item{category}{Show only results matching this category}

\item{from_date}{Start time window for search. Date format is YYYY-MM-DD}

\item{to_date}{Ending time window for search. Date format is YYYY-MM-DD}

\item{mine}{Browse only articles owned by user. default is FALSE. Not functional. Use \code{\link{fs_browse}} instead.}

\item{public_only}{(for use with mine=TRUE only) browse only my public articles. default is FALSE}

\item{private_only}{(for use with mine=TRUE only) browse only my private articles. default is FALSE}

\item{drafts_only}{(for use with mine=TRUE only) browse only my draft articles. default is FALSE}

\item{session}{(optional) the authentication credentials from \code{\link{fs_auth}}. If not provided, will attempt to load from cache as long as figshare_auth has been run.}

\item{base}{the API access url}

\item{debug}{logical, enable debugging mode}
}
\value{
output of PUT request (invisibly)
}
\description{
Search function that will filter on matching timestamp, author, title, description, tag, category, and date range.  Query searches against matches in any metadata field.  Full-text searches coming soon.
}
\examples{
\dontrun{
fs_search("Boettiger") 
fs_search("Boettiger", author = "Carl")
fs_search("Boettiger", author = "Carl", from="2014-01-01")
fs_search("Boettiger", author = "Carl", from="2014-01-01",
          category = "Evolutionary Biology")

} 
}
\references{
\url{http://api.figshare.com/docs/howto.html#q-search}
}
\seealso{
\\code{\link{fs_auth}} \code{\link{fs_browse}}
}
\author{
Carl Boettiger \email{cboettig@gmail.com}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fs_auth.R
\name{fs_auth}
\alias{fs_auth}
\title{Figshare authentication via OAuth1.0 using httr}
\usage{
fs_auth(cKey = getOption("FigshareKey", NULL),
  cSecret = getOption("FigsharePrivateKey", NULL),
  token = getOption("FigshareToken", NULL),
  token_secret = getOption("FigsharePrivateToken", NULL))
}
\arguments{
\item{cKey}{optional argument for the consumer key.  See details.}

\item{cSecret}{optional argument for the consumer secret key. See details.}

\item{token}{optional argument for the user-specific token.  See details.}

\item{token_secret}{Optional argument to provide a secret token assigned
to the user, rather than let fs_auth() automatically handle authentication. 
See details.}
}
\value{
OAuth credential (invisibly).  The credential is also written to the enivronment "FigshareAuthCache", which is created when the package is loaded.  All functions needing authentication can then access the credential without it being explicitly passed in the function call. If authentication fails, returns the failing GET response for debugging.
}
\description{
Figshare authentication via OAuth1.0 using httr
}
\details{
Explicit calls to fs_auth() are usually not needed,
as the function is called automatically by all other functions that
need authentication.  As of version 0.3, no arguments are needed as
authentication is done online, and fs_auth() will not attempt to load
keys stored in options.  

By default, the function will use the application's consumer key and
consumer secret key, rather than expecting the user to create their own
application.  The user-specific tokens will then be generated and locally
cached for use between sessions, if indicated by the interactive options. 
For details, see httr oauth1.0_token documentation.  

If for some reason a user would rather provide there token and secret token
as before this is still supported using the same arguments.  Users wanting to
have their own app can provide cKey and cSecret arguments too, but this
is provided primarily for backwards compatibility with older versions.  It
is expected that most users will leave the keys as NULL.
}
\examples{
\dontrun{
fs_auth()
}
}
\references{
\url{http://api.figshare.com}
}
\author{
Carl Boettiger \email{cboettig@gmail.com}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fs_embed.R
\name{fs_embed}
\alias{fs_embed}
\title{Upload a figure to figshare and return the url}
\usage{
fs_embed(file)
}
\arguments{
\item{file}{path to an image file}
}
\value{
a url to the image file
}
\description{
Upload a figure to figshare and return the url
}
\details{
use with opts_knit$set(upload.fn = fs_embed)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/summary_fs_details.R
\name{summary_fs_details}
\alias{summary_fs_details}
\title{Collect metadata from details call}
\usage{
summary_fs_details(fs_details_obj)
}
\arguments{
\item{fs_details_obj}{object}
}
\description{
Collect metadata from details call
}
\examples{
\dontrun{
fs_auth()
my_article <- fs_details("138")
summary_fs_details(my_article)
}
}
\references{
\url{http://api.figshare.com}
}
\author{
Edmund Hart \email{edmund.m.hart@gmail.com}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fs_upload.R
\name{fs_upload}
\alias{fs_upload}
\title{Upload file to an article}
\usage{
fs_upload(article_id, file, session = fs_get_auth())
}
\arguments{
\item{article_id}{an article id number or a character string (or list) of numbers}

\item{file}{path to file to upload, or character string (or list) of files (paths)}

\item{session}{the authentication credentials from \code{\link{fs_auth}} (optional)}
}
\description{
Upload file to an article
}
\details{
Articles may be draft, private or public but all uploads are saved as draft changes - the canonical public version of the deposit is not updated. To update the public version of the repository, use \code{\link{fs_make_public}}. Only articles of type "fileset" can have multiple files uploaded.

If only a single id number is given but a character string of files is given,
then be sure that the id corresponds to an object of type "fileset".  If article_id list
has more than one id, then there must be a corresponding file path for each id.
}
\examples{
\dontrun{
id <- fs_create("Title", "description", "figure")
fs_upload(id, "file.png")
} 
}
\references{
\url{http://api.figshare.com}
}
\seealso{
\code{\link{fs_auth}}
}
\author{
Carl Boettiger \email{cboettig@gmail.com}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fs_add_categories.R
\name{fs_add_categories}
\alias{fs_add_categories}
\title{Add a category to article}
\usage{
fs_add_categories(article_id, category_id, session = fs_get_auth(),
  debug = FALSE)
}
\arguments{
\item{article_id}{the id number of the article}

\item{category_id}{is a vector of integers corresponding to categories or a vector of category names}

\item{session}{(optional) the authentication credentials from \code{\link{fs_auth}}. If not provided, will attempt to load from cache as long as figshare_auth has been run.}

\item{debug}{return PUT results visibly?}
}
\value{
output of PUT request (invisibly)
}
\description{
Add a category to article
}
\examples{
\dontrun{
fs_add_categories(138, "Ecology")
}
}
\references{
\url{http://api.figshare.com}
}
\seealso{
\code{\link{fs_auth}}
}
\author{
Edmund Hart \email{edmund.m.hart@gmail.com}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fs_delete.R
\name{fs_delete}
\alias{fs_delete}
\title{Delete article (private or drafts only) or attached file}
\usage{
fs_delete(article_id, file_id = NULL, session = fs_get_auth(),
  debug = FALSE)
}
\arguments{
\item{article_id}{the id number of the article}

\item{file_id}{the id number of the file, if removing an attached file from a fileset.
file_id defaults to NULL, removing the entire draft or private article.}

\item{session}{(optional) the authentication credentials from \code{\link{fs_auth}}. If not provided, will attempt to load from cache as long as figshare_auth has been run.}

\item{debug}{display return value of request?}
}
\value{
output of DELETE request (invisibly)
}
\description{
Delete article (private or drafts only) or attached file
}
\examples{
\dontrun{
fs_delete(123)

## Delete all attachments in the second-most-recent entry in my library
my_lib <- fs_browse(mine=TRUE)
article_id <- my_lib[[2]]$article_id
file_ids <- sapply(my_lib[[2]]$files, `[[`, "id")
sapply(file_ids, function(id) fs_delete(article_id, id))
}
}
\references{
\url{http://api.figshare.com}
}
\seealso{
\code{\link{fs_auth}}
}
\author{
Carl Boettiger \email{cboettig@gmail.com}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fs_update.R
\name{fs_update}
\alias{fs_update}
\title{Update article title, description, or type}
\usage{
fs_update(article_id, title = NA, description = NA, type = NA,
  mine = TRUE, session = fs_get_auth(), debug = FALSE)
}
\arguments{
\item{article_id}{the id number of the article}

\item{title}{for the article (to replace original title)}

\item{description}{of the article (replaces original designation)}

\item{type}{one of: dataset, figure, media, poster, or paper (replaces original designation)}

\item{mine}{Set to \code{TRUE} if it refers to an item on your own account}

\item{session}{(optional) the authentication credentials from \code{\link{fs_auth}}.}

\item{debug}{return httr PUT request visibly?}
}
\value{
output of PUT request (invisibly)
}
\description{
Updates the article title, description or type. If any is not specified, it will remain unchanged.
}
\details{
Updates the title, description, and type of an article.
}
\examples{
\dontrun{
fs_update(138, title = "New title") 
}
}
\references{
\url{http://api.figshare.com}
}
\seealso{
\code{\link{fs_auth}}, \code{\link{fs_add_tags}}
}
\author{
Carl Boettiger \email{cboettig@gmail.com}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fs_add_authors.R
\name{fs_author_ids}
\alias{fs_author_ids}
\title{Get Author IDs from names}
\usage{
fs_author_ids(authors, session = fs_get_auth(), graphics = FALSE)
}
\arguments{
\item{authors}{a list/vector of authors (not a character string)}

\item{session}{(optional) the authentication credentials from \code{\link{fs_auth}}. If not provided, will attempt to load from cache.}

\item{graphics}{logical (default False) use graphic input to disambiguate?}
}
\value{
a list of author id numbers, or NULLS where ids cannot be found.
}
\description{
Take an author list, search for each author and return their FigShare ID.  
If no author is found, call fs_create_author and return the ID.
If multiple matches are found, allow user to choose interactively
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fs_details.R
\name{fs_details}
\alias{fs_details}
\title{Get details for an article}
\usage{
fs_details(article_id, mine = is_mine(article_id),
  session = fs_get_auth(), show_versions = FALSE, version = NULL,
  debug = FALSE)
}
\arguments{
\item{article_id}{number}

\item{mine}{logical (default FALSE). Set to true to see article details for your own non-public articles}

\item{session}{the authentication credentials from \code{\link{fs_auth}}}

\item{show_versions}{logical, show what versions are available}

\item{version}{show a given version number}

\item{debug}{logical, enable debugging mode?}
}
\description{
Get details for an article
}
\examples{
\dontrun{
fs_details(138)
}
}
\references{
\url{http://api.figshare.com}
}
\seealso{
\code{\link{fs_auth}}
}
\author{
Carl Boettiger \email{cboettig@gmail.com}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fs_browse.R
\name{fs_browse}
\alias{fs_browse}
\title{Browse articles}
\usage{
fs_browse(mine = TRUE, public_only = FALSE, private_only = FALSE,
  drafts_only = FALSE, session = fs_get_auth(),
  base = "http://api.figshare.com/v1", query = NA, debug = FALSE)
}
\arguments{
\item{mine}{Logical, show only my (authenticated user's) articles. Defaults to TRUE.}

\item{public_only}{(for use with mine=TRUE only) browse only my public articles. default is FALSE}

\item{private_only}{(for use with mine=TRUE only) browse only my private articles. default is FALSE}

\item{drafts_only}{(for use with mine=TRUE only) browse only my draft articles. default is FALSE}

\item{session}{(optional) the authentication credentials from \code{\link{fs_auth}}. If not provided, will attempt to load from cache as long as figshare_auth has been run.}

\item{base}{the API access url}

\item{query}{a search query term (equivalent to calling fs_search)}

\item{debug}{enable debugging mode}
}
\value{
output of PUT request (invisibly)
}
\description{
Browse can be set to all public articles, the users own articles, 
Browse can filter on matching timestamp, author, title, description, tag, category, and date range.
}
\examples{
\dontrun{
fs_browse() 
}
}
\references{
\url{http://api.figshare.com/docs/howto.html#q-search}
}
\seealso{
\code{\link{fs_auth}}
}
\author{
Carl Boettiger \email{cboettig@gmail.com}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fs_category_list.R
\name{fs_category_list}
\alias{fs_category_list}
\title{List all categories}
\usage{
fs_category_list(debug = FALSE)
}
\arguments{
\item{debug}{enable debugging}
}
\value{
a table of all the categories
}
\description{
List all categories
}
\examples{
\dontrun{
fs_categories_list()
}
}
\references{
\url{http://api.figshare.com}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fs_download.R
\name{fs_download}
\alias{fs_download}
\title{Get details for an article}
\usage{
fs_download(article_id, urls_only = TRUE, mine = is_mine(article_id),
  session = fs_get_auth(), show_versions = FALSE, version = NULL,
  ...)
}
\arguments{
\item{article_id}{number}

\item{urls_only}{logical (default TRUE) to only return the URLs to the 
downloadable objects but do not call download.file.  If FALSE, will download files}

\item{mine}{logical (default FALSE). Set to true to see article details for your own non-public articles}

\item{session}{the authentication credentials from \code{\link{fs_auth}}}

\item{show_versions}{logical, show what versions are available}

\item{version}{show a given version number}

\item{...}{additional arguments to \code{\link{download.file}}}
}
\description{
Get details for an article
}
\examples{
\dontrun{
url <- fs_download(90818)
data <- read.csv(url)
articles <- fs_search("SciFund")
ids <- fs_ids(articles)
fs_download(ids, urls_only=FALSE)
}
}
\references{
\url{http://api.figshare.com} \url{https://github.com/ropensci/rfigshare}
}
\seealso{
\code{\link{fs_auth}} \code{\link{download.file}}
}
\author{
Carl Boettiger \email{cboettig@gmail.com}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fs_new_article.R
\name{fs_new_article}
\alias{fs_new_article}
\title{Create a FigShare article.}
\usage{
fs_new_article(title, description, type = c("dataset", "figure", "media",
  "poster", "paper", "fileset"), authors = NA, categories = NA,
  tags = NA, links = NA, files = NA, visibility = c("draft",
  "private", "public"), session = fs_get_auth())
}
\arguments{
\item{title}{for the article, see \code{\link{fs_create}} for details.}

\item{description}{of the article, see \code{\link{fs_create}} for details.}

\item{type}{one of: dataset, figure, media, poster, or paper, see \code{\link{fs_create}} for details.}

\item{authors}{Orded list of authors for the article, see \code{\link{fs_add_authors}} for details}

\item{categories}{list of categories or category id numbers, see \code{\link{fs_add_categories}} for details.}

\item{tags}{list of tags, see \code{\link{fs_add_tags}} for details.}

\item{links}{list of links to add, see \code{\link{fs_add_links}} for details}

\item{files}{path to the files to add, see \code{\link{fs_upload}} for details}

\item{visibility}{one of "draft", "private" or "public".  A draft document can still be edited and modified.
A public document is visible to everyone and cannot be deleted (though additional authors to the work can still "claim" their authorship).}

\item{session}{(optional) credentials, see \code{link{fs_auth}}}
}
\value{
article id
}
\description{
fs_new_article is a wrapper around many other rfigshare commands to provide convenient posting.
}
\examples{
\dontrun{
write.csv(mtcars, "mtcars.csv")
id <- fs_new_article(title="A Test of rfigshare", 
                    description="This is a test of the fs_new_article function and related 
                    methods", 
                    type="dataset", 
                    authors=c("Karthik Ram", "Scott Chamberlain"), 
                    tags=c("ecology", "openscience"), 
                    categories="Ecology", 
                    links="http://ropensci.org", 
                    files="mtcars.csv",
                    visibility="private")
}
}
\references{
\url{http://api.figshare.com}
}
\seealso{
\code{\link{fs_auth}}, \code{\link{fs_add_categories}}, \code{\link{fs_add_authors}}, \code{\link{fs_add_tags}}, \code{\link{fs_add_links}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fs_make_private.R
\name{fs_make_private}
\alias{fs_make_private}
\title{Make an article private (draft only?)}
\usage{
fs_make_private(article_id, session = fs_get_auth())
}
\arguments{
\item{article_id}{the id number of the article}

\item{session}{(optional) the authentication credentials from \code{\link{fs_auth}}. If not provided, will attempt to load from cache as long as \code{\link{fs_auth}} has been run.}
}
\value{
output of PUT request (invisibly)
}
\description{
Make an article private (draft only?)
}
\examples{
\dontrun{
fs_make_private(123)
}
}
\references{
\url{http://api.figshare.com}
}
\seealso{
\code{\link{fs_auth}}
}
\author{
Carl Boettiger \email{cboettig@gmail.com}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fs_create.R
\name{fs_create}
\alias{fs_create}
\title{Create a FigShare article (draft)}
\usage{
fs_create(title, description, type = c("dataset", "figure", "media",
  "poster", "paper", "fileset"), session = fs_get_auth(),
  debug = FALSE)
}
\arguments{
\item{title}{for the article}

\item{description}{of the article}

\item{type}{one of: dataset, figure, media, poster, paper or fileset. (Only filesets can have multiple uploaded files attached).}

\item{session}{the authentication credentials from \code{\link{fs_auth}}}

\item{debug}{print full post call return}
}
\value{
article id
}
\description{
Articles must be created with \code{\link{fs_create}}
with essential metadata.  Then you can add files with
\code{\link{fs_upload}}, add categories, tags or authors
with \code{\link{fs_add_categories}} or \code{\link{fs_add_tags}}
\code{\link{fs_add_authors}}.  Authors not registered with a FigShare
id can be created with \code{\link{fs_create_author}}.  You can
edit the original metadata with \code{\link{fs_update}}.
Finally, release the article as either private or public with
\code{\link{fs_make_private}} or \code{\link{fs_make_public}}.
Before creating the article, you must authenticate using
\code{\link{fs_auth}}.
}
\examples{
\dontrun{
fs_create("My Title", "A description of the object", "dataset")
}
}
\references{
\url{http://api.figshare.com}
}
\seealso{
\code{\link{fs_auth}}, \code{\link{fs_upload}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fs_add_links.R
\name{fs_add_links}
\alias{fs_add_links}
\title{Add link to article}
\usage{
fs_add_links(article_id, link, session = fs_get_auth(), debug = FALSE)
}
\arguments{
\item{article_id}{the id number of the article}

\item{link}{the url you wish to add (can be list of urls)}

\item{session}{(optional) the authentication credentials from \code{\link{fs_auth}}. If not provided, will attempt to load from cache as long as authentication has been run.}

\item{debug}{logical, should function return details of PUT request?}
}
\value{
output of PUT request (invisibly)
}
\description{
Adds url links to the metadata of an article
}
\examples{
\dontrun{
fs_add_links(138, list("http://carlboettiger.info", "http://ropensci.org")) 
}
}
\references{
\url{http://api.figshare.com}
}
\seealso{
\code{\link{fs_auth}}
}
\author{
Carl Boettiger \email{cboettig@gmail.com}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_to_filename.R
\name{plot_to_filename}
\alias{plot_to_filename}
\title{Convienence function to save a ggplot2 plot, and return its filename.}
\usage{
plot_to_filename(plotobj, filename, path = ".")
}
\arguments{
\item{plotobj}{ggplot2 plot object (should add support for base plots too)}

\item{filename}{Filename, don't include the file type extension.}

\item{path}{Path where you want to save the file.}
}
\value{
A file name, to use in fs_upload
}
\description{
Convienence function to save a ggplot2 plot, and return its filename.
}
\examples{
\dontrun{ 
# include in your fs_upload call
library(ggplot2)
p <- qplot(mpg, wt, data=mtcars)
plott <- fs_create(title="my title", description="some description", type="figure")
fs_add_categories(plott, "Ecology")
fs_upload(plott, plot_to_filename(p, "myfilename", "~"))
}
}
\seealso{
\code{\link{fs_upload}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fs_author_search.R
\name{fs_author_search}
\alias{fs_author_search}
\title{Search for an author}
\usage{
fs_author_search(author, session = fs_get_auth(), debug = FALSE)
}
\arguments{
\item{author}{a string to search for (name, can include spaces)}

\item{session}{(optional) the authentication credentials from \code{\link{fs_auth}}. If not provided, will attempt to load from cache as long as figshare_auth has been run.}

\item{debug}{toggle debugging mode}
}
\value{
output of PUT request (invisibly)
}
\description{
Function to search for authors
}
\examples{
\dontrun{
fs_author_search("Boettiger") 
} 
}
\references{
\url{http://api.figshare.com}
}
\seealso{
\code{\link{fs_auth}}
}
\author{
Carl Boettiger \email{cboettig@gmail.com}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fs_add_tags.R
\name{fs_add_tags}
\alias{fs_add_tags}
\title{Add a tag to an article}
\usage{
fs_add_tags(article_id, tag, session = fs_get_auth(), debug = FALSE)
}
\arguments{
\item{article_id}{the id number of the article to create}

\item{tag}{name of the tag to add (or list of tags)}

\item{session}{the authentication credentials from \code{\link{fs_auth}}}
}
\value{
output of PUT request (invisibly)
}
\description{
Add a tag to an article
}
\examples{
\dontrun{
 fs_add_tag(138, "phylogenetics") 
}
}
\references{
\url{http://api.figshare.com}
}
\seealso{
\code{\link{fs_auth}}
}
\author{
Carl Boettiger \email{cboettig@gmail.com}
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fs_add_categories.R
\name{fs_cat_to_id}
\alias{fs_cat_to_id}
\title{Helper function that matches string categories to id's}
\usage{
fs_cat_to_id(category_id)
}
\arguments{
\item{category_id}{Must be a valid category string, regardless of case}
}
\value{
a vector of integers corresponding to valid figshare categories
}
\description{
Helper function that matches string categories to id's
}
\references{
\url{http://api.figshare.com}
}
\author{
Edmund Hart \email{edmund.m.hart@gmail.com}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fs_embed.R
\name{fs_image_url}
\alias{fs_image_url}
\title{get image url}
\usage{
fs_image_url(id, debug = FALSE)
}
\arguments{
\item{id}{a (public) figshare figure id number}

\item{debug}{logical, enable debugging mode?}
}
\value{
a url to the image file
}
\description{
get image url
}
\details{
this is currently an unstable hack of html parsing
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fs_add_authors.R
\name{fs_add_author}
\alias{fs_add_author}
\title{Add an author to an article by ID number}
\usage{
fs_add_author(article_id, author_id, session = fs_get_auth())
}
\arguments{
\item{article_id}{id number of an article on figshare}

\item{author_id}{the id number of a registered figshare user (see \code{\link{fs_author_search}})}

\item{session}{(optional) the authentication credentials from \code{\link{fs_auth}}. If not provided, will attempt to load from cache.}
}
\description{
Add an author to an article by ID number
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fs_add_authors.R
\name{fs_add_authors}
\alias{fs_add_authors}
\title{Add author to an article}
\usage{
fs_add_authors(article_id, authors, session = fs_get_auth(),
  create_missing = TRUE, debug = FALSE)
}
\arguments{
\item{article_id}{id number of an article on figshare}

\item{authors}{a list or character string of authors or author id numbers (or mixed).}

\item{session}{(optional) the authentication credentials from \code{\link{fs_auth}}. 
If not provided, will attempt to load from cache as long as figshare_auth has been run.}

\item{create_missing}{(logical) Attempt to create authors not already registered on FigShare?}

\item{debug}{return the httr result visibly?}
}
\value{
adds the requested authors to the given article
}
\description{
Add author to an article
}
\examples{
\dontrun{
 fs_add_authors("138", list("Scott Chamberlain", "Karthik Ram"))
 fs_add_authors("138", c("Scott Chamberlain", "Karthik Ram"))
 fs_add_authors("138", list("Scott Chamberlain", "97306"))
 fs_add_authors("138", list("Scott Chamberlain", 97306))
 fs_add_authors(138, 97306)
} 
}
\author{
Carl Boettiger \email{cboettig@gmail.com}
}
