[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1296448.svg)](https://doi.org/10.5281/zenodo.1296448)

# staphopia-paper
Code used for Staphopia manuscript
---
title: "Results Section: Public Metadata"
output:
  pdf_document: default
---
```{r}
library(staphopia)
library(dplyr)
library(ggplot2)
library(reshape2)
```

## Aggregating Data For Public Samples
First we'll get all publicly available *S. aureus* samples.

```{r}
ps <- get_public_samples()
```

We now have `r toString(nrow(ps))` samples to work with. Next we will acquire
metadata associated with each sample.

We will also get information pertaining to submissions by year and how any publication links were made.
```{r}
submissions <- get_submission_by_year()
publication_links <- get_publication_links()
```

Next we are going to pull down any metadata associated with the public samples.
```{r}
metrics <- merge(
    ps, 
    get_metadata(ps$sample_id),
    by='sample_id'
)
```

We are now going to add two columns `rank_name` and `year`.
```{r}
metrics$year <- sapply(
    metrics$first_public,
    function(x) {
        strsplit(x, "-")[[1]][1]
    }
)

metrics$rank_name <- ifelse(
    metrics$rank.x == 3,
    'Gold',
    ifelse(
        metrics$rank.x == 2,
        'Silver',
        'Bronze'
    )
)
```

### Publication Information

#### Summary
Here are details looking at total submissions and their publication status.
```{r}
t(submissions[submissions$year == max(submissions$year),])
```

Here is information on how publication links were made.
```{r}
t(publication_links)
```

There are 6 rows and their names are as follows:

1. elink: Number samples linked to a PubMed ID identified from eLink
2. text: Number samples linked to a PubMed ID identified from text mining (not through eLink)
3. elink_pmid: Number of PubMed IDs identified from eLink
4. text_pmid: Number of PubMed IDs identified from text mining (not through eLink)
5. total: Total number of samples associated with a PubMed ID
6. total_pmid: Total number of PubMed IDs associated with published samples

##### Percent of Samples Published
```{r}
stats <- submissions[submissions$year == max(submissions$year),]
stats$overall_published / stats$overall * 100
```

#### Published vs Unpublished Submissions Per Year
```{r fig.width=12, fig.asp=0.5}
melted <- melt(submissions, id=c('year'),
               measure.vars = c('published', 'unpublished'))
melted$title <- ifelse(melted$variable == 'published', 'Published', 'Unpublished')
p <- ggplot(data=melted, aes(x=year, y=value, fill=title)) +
    xlab("Year") +
    ylab("Count") +
    geom_bar(stat='identity', position='dodge') +
    geom_text(aes(label=value), vjust = -0.5, position = position_dodge(.9)) +
    scale_fill_manual(values=c("#2ca25f", "#5ab4ac")) +
    scale_x_continuous(breaks = round(
        seq(min(submissions$year), max(submissions$year), by = 1), 1
    )) +
    theme_bw() +
    theme(axis.text=element_text(size=12),
          axis.title=element_text(size=14,face="bold"),
          legend.title = element_blank())
    

p
```

#### Overall Published vs Unpublished Submissions
```{r Figure Pubs, fig.width=12, fig.asp=0.5}
melted <- melt(submissions, id=c('year'),
               measure.vars = c('overall_published', 'overall_unpublished'))
melted$title <- ifelse(melted$variable == 'overall_published', 'Published', 'Unpublished')
p <- ggplot(data=melted, aes(x=year, y=value, fill=title)) +
    xlab("Year") +
    ylab("Cumulative Count") +
    geom_bar(stat='identity', position='dodge') +
    geom_text(aes(label=value), vjust = -0.5, position = position_dodge(.9)) +
    scale_fill_manual(values=c("#2ca25f", "#5ab4ac")) +
    scale_x_continuous(breaks = round(
        seq(min(submissions$year), max(submissions$year), by = 1), 1
    )) +
    theme_bw() +
    theme(axis.text=element_text(size=12),
          axis.title=element_text(size=14,face="bold"),
          legend.title = element_blank())
    

p
```

#### Overall Published vs Unpublished Submissions
```{r Figure 8 Pubs, fig.width=12, fig.asp=0.5}
melted <- melt(submissions, id=c('year'),
               measure.vars = c('overall_published', 'overall_unpublished'))
melted$title <- ifelse(melted$variable == 'overall_published', 'Published', 'Unpublished')
melted$final <- ifelse(melted$year == 2017, melted$value, "") 

p <- ggplot(data=melted, aes(x=year, y=value, fill=title, label=final)) +
    xlab("Year") +
    ylab("Cumulative Count") +
    geom_bar(stat='identity', position='stack') +
    # geom_text(position = position_stack(vjust = 0.50)) +
    geom_text(aes(year, overall + 1000, label = overall, fill = NULL), data = submissions) +
    scale_fill_manual(values=c("#2ca25f", "#5ab4ac")) +
    scale_x_continuous(breaks = round(
        seq(min(submissions$year), max(submissions$year), by = 1), 1
    )) +
    theme_bw() +
    theme(axis.text=element_text(size=12),
          axis.title=element_text(size=14,face="bold"),
          legend.title = element_blank())
    

p

# Output plot to PDF and PNG
staphopia::write_plot(p, paste0(getwd(), '/../figures/figure-01-submission-published-per-year'))
```

### Metadata Information

#### Number of Samples With A Collection Date
```{r}
has_collection_date <- nrow(metrics[metrics$collection_date != "",])
paste0(has_collection_date," (", has_collection_date / nrow(metrics) * 100, " %)")
```

#### Number of Samples With A Location Information
```{r}
has_location <- nrow(metrics[metrics$location != "unknown/missing",])
paste0(has_location," (", has_location / nrow(metrics) * 100, " %)")
```

#### Number of Locations
```{r}
nrow(as.data.frame(table(metrics[metrics$location != "unknown/missing",]$location)))
```

#### Countries
```{r}
country_data <- as.data.frame(table(
    metrics[(metrics$country != "unknown/missing" ) & (metrics$country != ""),]$country
))
colnames(country_data) <- c("Country", "total")
country_data <- arrange(country_data, desc(total))
country_data
```

#### Number of Countries
```{r}
paste0(nrow(country_data), " countries, represented by ", sum(country_data$total), " samples")
```

#### Number of Samples With Isolation Source
```{r}
has_source <- nrow(metrics[metrics$isolation_source != "",])
paste0(has_source," (", has_source / nrow(metrics) * 100, " %)")
```

#### Isolation Sources
```{r}
df <- as.data.frame(table(substr(tolower(
    metrics[metrics$isolation_source != "",]$isolation_source), 1, 50
)))
df[order(-df$Freq),]
```

#### Number of Isolation Sources
```{r}
nrow(as.data.frame(table(tolower(
    metrics[metrics$isolation_source != "",]$isolation_source
))))
```

# Session Info
```{r}
sessionInfo()
```
---
title: "Results Section: Public Sequencing Metrics"
output:
  pdf_document: default
---
```{r}
library(staphopia)
library(ggplot2)
library(reshape2)
```

## Aggregating Data For Public Samples
First we'll get all publicly available *S. aureus* samples.

```{r}
ps <- get_public_samples()
```

We will also get information pertaining to submissions and ranks by year.
```{r}
submissions <- get_submission_by_year(all = TRUE)
ranks <- get_rank_by_year()
```

We now have `r toString(nrow(ps))` samples to work with. Next we will acquire
metadata, sequencing stats and assembly stats associated with each sample.

```{r}
metrics <- merge(
    ps,
    merge(
        get_assembly_stats(ps$sample_id),
        merge(
            get_metadata(ps$sample_id),
            get_sequence_quality(ps$sample_id, stage='cleanup'),
            by='sample_id'
        ),
        by='sample_id'
    ),
    by='sample_id'
)
```

We are now going to add two columns `rank_name` and `year`.
```{r}
metrics$year <- sapply(
    metrics$first_public,
    function(x) {
        strsplit(x, "-")[[1]][1]
    }
)

metrics$rank_name <- ifelse(
    metrics$rank.x == 3,
    'Gold',
    ifelse(
        metrics$rank.x == 2,
        'Silver',
        'Bronze'
    )
)
```

## Visualizing Metrics
The following sections will be plots to visualize relationships in the data.

### By Year Plots
#### Submissions Per Year
```{r fig.width=12, fig.asp=0.4}
p <- ggplot(data=submissions, aes(x=year, y=count)) +
    xlab("Year") +
    ylab("Count") +
    geom_bar(stat='identity') +
    geom_text(aes(label=count), vjust = -0.5) +
    scale_x_continuous(breaks = round(
        seq(min(submissions$year), max(submissions$year), by = 1),1
    )) +
    theme_bw() +
    theme(axis.text=element_text(size=12),
          axis.title=element_text(size=14,face="bold"))
p
```

#### Overall Submissions
```{r Figure 1, fig.asp=0.4, fig.width=12}
p <- ggplot(data=submissions, aes(x=year, y=overall)) +
    xlab("Year") +
    ylab("Cumulative Count") +
    geom_bar(stat='identity') +
    geom_text(aes(label=overall), vjust = -0.5) +
    scale_x_continuous(breaks = round(
        seq(min(submissions$year), max(submissions$year), by = 1), 1
    )) +
    theme_bw() +
    theme(axis.text=element_text(size=12),
          axis.title=element_text(size=14,face="bold"))
p
```

#### Submission Ranks
```{r Figure 4 BSG, fig.width=12, fig.asp=0.5}
melted <- melt(ranks, id=c('year'),
               measure.vars = c('bronze', 'silver', 'gold'))
melted$title <- ifelse(melted$variable == 'gold', 'Gold', 
                       ifelse(melted$variable == 'silver', 'Silver', 'Bronze'))
melted$rank <- ifelse(melted$variable == 'gold', 3, 
                      ifelse(melted$variable == 'silver', 2, 1))
p <- ggplot(data=melted, aes(x=year, y=value, fill=title, group=rank, label=title)) +
    xlab("Year") +
    ylab("Count") +
    geom_bar(stat='identity', position='dodge') +
    geom_text(aes(label=value), vjust = -0.5, position = position_dodge(.9)) +
    scale_fill_manual(values=c("#CD7F32", "#D4AF37", "#C0C0C0")) +
    scale_x_continuous(breaks = round(
        seq(min(ranks$year), max(ranks$year), by = 1), 1
    )) +
    theme_bw() +
    theme(axis.text=element_text(size=12),
          axis.title=element_text(size=14,face="bold"),
          legend.title = element_blank())
    

p

# Output plot to PDF and PNG
staphopia::write_plot(p, paste0(getwd(), '/../figures/figure-03-rank-per-year'))
```

#### Assembly Size
```{r fig.width=12, fig.asp=0.5}
p <- ggplot(metrics, aes(x = year, y = total_contig_length)) +
    geom_boxplot()
p
```

#### Total Contigs (smaller is better)
```{r fig.width=12, fig.asp=0.5}
p <- ggplot(metrics, aes(x = year, y = total_contig)) +
    geom_boxplot()
p
```

#### N50
```{r fig.width=12, fig.asp=0.5}
p <- ggplot(metrics, aes(x = year, y = n50_contig_length)) +
    geom_boxplot()
p
```

#### Mean Contig Length
```{r fig.width=12, fig.asp=0.5}
p <- ggplot(metrics, aes(x = year, y = mean_contig_length)) +
    geom_boxplot()
p
```

#### Max Contig Length
```{r fig.width=12, fig.asp=0.5}
p <- ggplot(metrics, aes(x = year, y = max_contig_length)) +
    geom_boxplot()
p
```

#### Mean Read Length
```{r fig.width=12, fig.asp=0.5}
p <- ggplot(metrics, aes(x = year, y = read_mean)) +
    geom_boxplot()
p
```

#### Mean Per-Read Quality Score
```{r fig.width=12, fig.asp=0.5}
p <- ggplot(metrics, aes(x = year, y = qual_mean)) +
    geom_boxplot()
p
```

#### Assembly Size Grouped By Rank
```{r fig.width=12, fig.asp=0.5}
p <- ggplot(metrics, aes(x = year, y = total_contig_length,
                         fill=rank_name, label=rank_name)) +
    geom_boxplot()
p
```

#### Total Contigs Grouped By Rank
```{r fig.width=12, fig.asp=0.5}
p <- ggplot(metrics, aes(x = year, y = total_contig,
                         fill=rank_name, label=rank_name)) +
    geom_boxplot()
p
```

#### N50 Grouped By Rank
```{r fig.width=12, fig.asp=0.5}
p <- ggplot(metrics, aes(x = year, y = n50_contig_length,
                         fill=rank_name, label=rank_name)) +
    geom_boxplot()
p
```

#### Mean Contig Length Grouped By Rank
```{r fig.width=12, fig.asp=0.5}
p <- ggplot(metrics, aes(x = year, y = mean_contig_length,
                         fill=rank_name, label=rank_name)) +
    geom_boxplot()
p
```

#### Max Contig Length Grouped By Rank
```{r fig.width=12, fig.asp=0.5}
p <- ggplot(metrics, aes(x = year, y = max_contig_length,
                         fill=rank_name, label=rank_name)) +
    geom_boxplot()
p
```

#### Mean Read Length Grouped By Rank
```{r fig.width=12, fig.asp=0.5}
p <- ggplot(metrics, aes(x = year, y = read_mean,
                         fill=rank_name, label=rank_name)) +
    geom_boxplot()
p
```

#### Mean Per-Read Quality Score Grouped By Rank
```{r fig.width=12, fig.asp=0.5}
p <- ggplot(metrics, aes(x = year, y = qual_mean,
                         fill=rank_name, label=rank_name)) +
    geom_boxplot()
p
```

### By Rank Plots
#### Assembly Size
```{r fig.width=12, fig.asp=0.5}
p <- ggplot(metrics, aes(x = rank_name, y = total_contig_length)) +
    geom_boxplot()
p
```

#### Total Contigs (smaller is better)
```{r fig.width=12, fig.asp=0.5}
p <- ggplot(metrics, aes(x = rank_name, y = total_contig)) +
    geom_boxplot()
p
```

#### N50
```{r fig.width=12, fig.asp=0.5}
p <- ggplot(metrics, aes(x = rank_name, y = n50_contig_length)) +
    geom_boxplot()
p
```

#### Mean Contig Length
```{r fig.width=12, fig.asp=0.5}
p <- ggplot(metrics, aes(x = rank_name, y = mean_contig_length)) +
    geom_boxplot()
p
```

#### Max Contig Length
```{r fig.width=12, fig.asp=0.5}
p <- ggplot(metrics, aes(x = rank_name, y = max_contig_length)) +
    geom_boxplot()
p
```

#### Mean Read Length
```{r fig.width=12, fig.asp=0.5}
p <- ggplot(metrics, aes(x = rank_name, y = read_mean)) +
    geom_boxplot()
p
```

#### Mean Per-Read Quality Score
```{r fig.width=12, fig.asp=0.5}
p <- ggplot(metrics, aes(x = rank_name, y = qual_mean)) +
    geom_boxplot()
p
```

#### Coverage
```{r fig.width=12, fig.asp=0.5}
p <- ggplot(metrics, aes(x = rank_name, y = coverage)) +
    geom_boxplot()
p
```

## Bronze Data
#### Coverage By Quality
```{r fig.width=12, fig.asp=0.5}
p <- ggplot(metrics[metrics$rank.x == 1,], aes(x = coverage, y = qual_mean)) +
    geom_point()
p
```

#### Coverage By Read Length
```{r fig.width=12, fig.asp=0.5}
p <- ggplot(metrics[metrics$rank.x == 1,], aes(x = coverage, y = read_mean)) +
    geom_point()
p
```


## Session Info
```{r}
sessionInfo()
```
---
title: "Results Section: Public Genetic Diversity"
output:
  pdf_document: default
---
```{r}
library(staphopia)
library(ggplot2)
library(reshape2)
library(scales)
```

In this section we will look into genetic diversity that has been sequenced in *Staphylococcus aureus*. In order to do so, we'll use variant counts, cgMLST and MLST as measures of diversity.


## Aggregating Data For Public Samples
First we'll get all publicly available *S. aureus* samples.

```{r}
ps <- get_public_samples()
```

## MLST
Next we will will use the MLST information has a measure of genitic diversity. In this case we are interested in the total number of unique sequence types sequenced. We'll use *get_st_by_year()* to get some basic stats about how many STs have been sequenced. We will also use *get_top_sequence_types()* to get each ST represented in the database and the total number of samples with each ST. (*Note: 5000 is just an arbitrarly large number to retreive all STs*)
```{r}
sequence_types <- get_st_by_year()
top_st <- get_top_sequence_types(5000)
colnames(sequence_types)
```

This gives us 38 columns for each year. These columns are:

1. year: The year.
2. unique: The Number of unique STs for a given year.
3. novel: Number of STs not sequenced previously.
4. assigned: Samples which a ST was determined.
5. assigned_agree: Samples in which each program that called an ST agreed in ST.
6. assigned_disagree: Samples in which programs did not each call the same ST.
7. unassigned: Samples which a ST was not determined.
8. unassigned_agree: Each program was unable to assign an ST.
9. unassigned_disagree: Samples in which no ST was determined, but each program does not agree
10. predicted_novel: Samples with a match to each Loci, but allele pattern does not exist.
11. all: Samples with an ST determined with agreement between each program.
12. partial: Samples with an ST determined with agreement between two programs.
13: ariba_blast: Samples with an ST determined with agreement between Ariba and BLAST.
14. mentalist_blast: Samples with an ST determined with agreement between MentaLiST and BLAST.
15. mentalist_ariba: Samples with an ST determined with agreement between MentaLiST and Ariba.
16. single: Samples with an ST determined by only a single program.
17. ariba: Samples with an ST determined by only Ariba.
18. mentalist: Samples with an ST determined by only MentaLiST.
19. blast: Samples with an ST determined by only BLAST.
20. count: Total number of samples in a given year.
21-38: overall_*X*: The cumulative totals of previous years for column *x*

#### Compare MLST Predictions
```{r}
mlst <- get_sequence_type(ps$sample_id)
metadata <- merge(
    ps, 
    get_metadata(ps$sample_id),
    by='sample_id'
)
metadata$year <- sapply(
    metadata$first_public,
    function(x) {
        strsplit(x, "-")[[1]][1]
    }
)

metadata$rank_name <- ifelse(
    metadata$rank == 3,
    'Gold',
    ifelse(
        metadata$rank == 2,
        'Silver',
        'Bronze'
    )
)
```


```{r}
mlst_temp <- merge(mlst, metadata[,c('sample_id', 'is_paired')], by='sample_id')
mlst_temp$is_paired <- ifelse(mlst_temp$is_paired == "", FALSE, TRUE)


mlst_temp$agreement <- paste0(
    ifelse(mlst$st == 0 | mlst$st == 0, '000',
        ifelse(mlst$mentalist == mlst$ariba & mlst$mentalist == mlst$blast, '111', 
            ifelse(mlst$mentalist == mlst$ariba, '110', 
                ifelse(mlst$mentalist == mlst$blast, '101',
                    ifelse(mlst$blast == mlst$ariba, '011',
                        ifelse(mlst$mentalist > 0, '100', 
                            ifelse(mlst$ariba > 0, '010',
                                ifelse(mlst$blast > 0, '001', '000')
                            )
                        )
                    )
                )
            )
        )
    )
)

# mentalist
# ariba
# blast
mlst_temp$agreement<- ifelse(mlst_temp$is_paired == TRUE, mlst_temp$agreement,
                             paste0(substr(mlst_temp$agreement, 1, 1), '-',
                                    substr(mlst_temp$agreement, 3, 3)))
table(mlst_temp$agreement)
platform <- metadata[,c('sample_id', 'instrument_model', 'study_accession', 'year', 'rank_name')]
mlst_temp <- merge(mlst_temp, platform, by='sample_id')
```

```{r}
table(mlst_temp[mlst_temp$agreement == '101',]$rank_name)
table(mlst_temp[mlst_temp$agreement == '011',]$rank_name)
table(mlst_temp[mlst_temp$agreement == '110',]$rank_name)
```

### PubMLST ST Counts
```{r}
st_counts <- merge(
    read.table('../data/pubmlst-counts.txt', header=TRUE, sep="\t"),
    staphopia <- top_st[top_st$st > 0,c('st', 'count')],
    by='st', all=TRUE
)
st_counts[is.na(st_counts)] <- 0
st_counts <- st_counts[st_counts$pubmlst_count > 0,]
nrow(st_counts)
nrow(st_counts[st_counts$count == 0,])
nrow(st_counts[st_counts$count >= 1,])
nrow(st_counts[st_counts$pubmlst_count == 1,])
nrow(st_counts[st_counts$pubmlst_count <= 2,])
table(st_counts[st_counts$count == 0,]$pubmlst_count)
summary(st_counts[st_counts$count > 0,]$pubmlst_count)
st_counts[st_counts$pubmlst_count > 10 & st_counts$count == 0,]
```

### Summary of MLST Diversity
#### Assignment Breakdown
```{r}
t(sequence_types[sequence_types$year == max(sequence_types$year),21:38])
```

#### Top STs
```{r}
top_st[1:10,]
```

This gives us 4 columns for each ST, in descending order based on the *count* column. In other words the most represented STs are seen first. These columns are:

1. st: The sequence type.
2. count: The number of samples with given ST.
3. percent: The percent of samples represented by given ST.
4. overall: The percent of samples represented by given ST and previous STs.

##### How many unique STs represented?
```{r}
nrow(top_st[top_st$st > 0,])
```

##### How many STs represented by a single sample?
```{r}
nrow(top_st[top_st$count == 1, ])
```


###  Visualizing MLST Diversity
The following sections will be plots to visualize relationships in the data.

#### Unique Sequence Types By Year
```{r fig.width=12, fig.asp=0.4}
p <- ggplot(data=sequence_types, aes(x=year, y=unique)) +
    xlab("Year") +
    ylab("Count") +
    geom_bar(stat='identity') +
    geom_text(aes(label=unique), vjust = -0.5) +
    scale_x_continuous(breaks = round(seq(min(sequence_types$year), 
                                          max(sequence_types$year), 
                                          by = 1),1)) +
    theme_bw() +
    theme(axis.text=element_text(size=12),
          axis.title=element_text(size=14,face="bold"))
p
```

#### Novel Sequence Types By Year
```{r fig.width=12, fig.asp=0.4}
p <- ggplot(data=sequence_types, aes(x=year, y=novel)) +
    xlab("Year") +
    ylab("Count") +
    geom_bar(stat='identity') +
    geom_text(aes(label=novel), vjust = -0.5) +
    scale_x_continuous(breaks = round(seq(min(sequence_types$year), 
                                          max(sequence_types$year), 
                                          by = 1),1)) +
    theme_bw() +
    theme(axis.text=element_text(size=12),
          axis.title=element_text(size=14,face="bold"))
p
```

#### Overall Novel Sequence Types By Year
```{r fig.width=12, fig.asp=0.4}
p <- ggplot(data=sequence_types, aes(x=year, y=overall_novel)) +
    xlab("Year") +
    ylab("Cumulative Count") +
    geom_bar(stat='identity') +
    geom_text(aes(label=overall_novel), vjust = -0.5) +
    scale_x_continuous(breaks = round(seq(min(sequence_types$year), 
                                          max(sequence_types$year), 
                                          by = 1),1)) +
    theme_bw() +
    theme(axis.text=element_text(size=12),
          axis.title=element_text(size=14,face="bold"))
p
```


#### Top 10 Sequence Types
```{r Figure 5 TopST, fig.width=12, fig.asp=0.45}
p <- ggplot(data=top_st[1:10,], aes(x=reorder(st, -count), y=count)) +
    xlab("Sequence Type") +
    ylab("Count") +
    geom_bar(stat="identity") +
    geom_text(aes(label=count), vjust = -0.5) + 
    theme_bw() +
    theme(axis.text=element_text(size=12),
          axis.title=element_text(size=14,face="bold"))
p
```

#### Total Allele Matches For Unassigned Samples
```{r fig.width=12, fig.asp=0.45}
allele_matches <- get_mlst_allele_matches(ps[ps$st == 0,]$sample_id)
df <- as.data.frame(table(allele_matches[allele_matches$matches < 7,]$matches))
colnames(df) <- c("matches", "count")

p <- ggplot(data=df, aes(x=matches, y=count)) +
    xlab("Matched Alleles") +
    ylab("Count") +
    geom_bar(stat="identity") +
    geom_text(aes(label=count), vjust = -0.5) + 
    theme_bw() +
    theme(axis.text=element_text(size=12),
          axis.title=element_text(size=14,face="bold"))
p
```

## cgMLST Patterns
Finally, we'll look at cgMLST as a measure of genetic diversity. We will use the *get_cgmlst()* function to get the cgMLST results for each Sample. This function might take a little while to retrieve all teh results.
```{r}
# USE_DEV to prevent timeout here until problem resolved
USE_DEV = TRUE

cgmlst <- get_public_cgmlst_patterns()
cgmlst$percent <- cgmlst$count / sum(cgmlst$total_samples)
cgmlst
```

This gives us two columns:

1. samples_in_pattern: The number of samples with a given cgMLST pattern.
2. count: The number patterns with a given number of samples.
3. total_samples: Number of samples represented by a row (samples_in_pattern * count)
4. percent: Percent of samples represented

For example, if samples_in_pattern is 100 and the count is 2. That means there are **2** (count=2) cgMLST patterns that are shared by **100 samples** (samples_in_count=100) each, representing a total of **200 samples** (count * samples_in_count).

### Total Number of Distinct cgMLST Patterns
```{r}
sum(cgmlst$count)
```

### How many shared cgMLST patterns?
```{r}
sum(cgmlst[cgmlst$samples_in_pattern > 1, ]$count)
```

### How many samples share a cgMLST pattern?
```{r}
sum(cgmlst[cgmlst$samples_in_pattern > 1, ]$total_samples)
```

### How many samples have a unique cgMLST pattern?
```{r}
cgmlst$percent <- cgmlst$count / sum(cgmlst$total_samples)
cgmlst[cgmlst$samples_in_pattern == 1, ]
```


## Session Info
```{r}
sessionInfo()
```
---
title: "Results Section: Supplementary subsampling"
output:
  pdf_document: default
---

## Read In The Data
```{r}
library(ggplot2)
library(dplyr)
options(scipen=999)

make_plot <- function (df, ylab) {
    p <- ggplot(data=df, aes(x=x, y=y)) +
        ylab(ylab) +
        xlab("Coverage") +
        geom_point(aes(color=color)) +    
        scale_x_continuous(breaks = seq(min(df$x), max(df$x), by = 20)) +
        theme_bw() +
        theme(axis.text=element_text(size=12),
              axis.title=element_text(size=14,face="bold"))
    return(p)
}

results <- read.table(
    "../data/supplementary-subsample/subsample-summary.txt", 
    header = TRUE,
    sep = "\t"
)
results <- results[results$mutations < 10000 & results$simulation != 'EF',]
colnames(results)
```

## Overview
We wanted to determine if at what level of coverage diminishing returns were observed in our analysis. This is important because high coverage sequences require more compuational resources (mostly in the form of memory) and take longer to process. Because we used Cancer Genomics Cloud (http://www.cancergenomicscloud.org/) to process each project, this was also important to reduce overall costs. We used runtime, assembly metrics and singleton kmer counts to determine a cutoff for coverage at which further coverage did not improve the results.


## Simulating Sequencing
We simulated HiSeq and MiSeq sequencing of the *S. aureus* N315 (NC_00274) reference genome with ART (Huang, W., Li, L., Myers, J.R., Marth, G.T., 2012. ART: a next-generation sequencing read simulator. Bioinformatics 28, 593–594.). Multiple coverages were simulated (see below) and processed through the Staphopia analysis pipeline on CGC.

We simulated multiple coverages:
```{r}
sort(unique(results$coverage))
```

## How does coverage affect run time?
```{r echo=FALSE, fig.asp=0.4, fig.width=12}
p <- make_plot(
    data.frame(x=results$coverage, y=results$runtime, color=results$simulation),
    "Total Run Time In Minutes"
)
p
```

In the plot above, there is evidence that increasing coverage leads to longer runtimes. Based on simulations MiSeq (MS) seqeuncing tended to take longer to process than HiSeq (HS). 

## How does coverage affect costs?
```{r fig.width=12, fig.asp=0.4}
p <- make_plot(
    data.frame(x=results$coverage, y=results$price, color=results$simulation),
    "Total Price For Analysis In USD"
)
p
```

The job cost on CGC is dependent on overall runtime. In the plot above, there is not much of a difference between 75x and 125x (~$0.125), but at 300x (~$0.175) it is about a $0.05 difference. For 44,000 genomes, it costs ~$5,500 to process genomes at a 75-125x coverage cutoff and ~$7,700 to process genomes at a 300x coverage cutoff. This is roughly a $2,200 difference in price.

## How does coverage affect assembly?
### Total number of contigs
```{r fig.width=12, fig.asp=0.4}
p <- make_plot(
    data.frame(x=results$coverage, y=results$total_contig, color=results$simulation),
    "Total Number of Assembled Contigs"
)
p
```


In the plot above, at 20x and 50x coverages there are more contigs that at >75x coverage, suggesting these coverages may not produce the best assembly. At 75x coverage and onwards, the total number of contigs does not change much. At >200x, there looks like a slight increase in total contigs.


### Total number of contigs greater than 200bp
```{r fig.width=12, fig.asp=0.4}
p <- make_plot(
    data.frame(x=results$coverage, y=results$total_contig_200bp, color=results$simulation),
    "Total Contigs >200bp Length"
)
p
```

Looking at the total number of contigs greater than 200bp, a similar pattern is oberseved. Again 20x and 50x may not produce the best assembly, and 75x coverage an onwards produce similar numbers.


### N50 contig length
```{r fig.width=12, fig.asp=0.4}
p <- make_plot(
    data.frame(x=results$coverage, y=results$n50_contig_length, color=results$simulation),
    "Assembly N50 Length"
)
p
```

Again, similar to the above to plots. THe exception is the N50 appears to level off around 100x instead of 75x.


## How does coverage affect total number of singleton kmers?
```{r fig.width=12, fig.asp=0.4}
p <- make_plot(
    data.frame(x=results$coverage, y=results$total_singleton, color=results$simulation),
    "Total Singleton k-mer Count"
)
p
```

In the plot above, it appears that further sequencing depth increases the number of observed singleton kmers. While there may be true singletons in the sequencing, many of these can be assumed to be due to sequencing errors. This suggests that at higher seqeuncing depths, there is greater need to correct erroneous reads that can affect analysis results. The effect is much greater in MiSeq than HiSeq, most likely due to the difference in error profiles used by ART.


## Conclusions
Overall using the metrics described above it appears coverage cutoff can be used without affecting the results of an analysis. Although, samples with 20x and 50x coverage had the lowest runtimes they produced assemblies that could be improved by further coverage. For samples with 150x or greater coverage produced assemblies similar to 75x-125x coverage, but took longer to process (increased costs) and also had more singleton kmers. This leaves 75x, 100x, and 125x coverages that were each very similar in runtime, costs, assemblies and kmers. Out of these three coverages, we arbitrary selected 100x as the default coverage cutoff for Staphopia analysis. 
---
title: "Results Section: Antibiotic Resistance Patterns"
output:
  pdf_document: default
---
```{r}
library(staphopia)
library(ggplot2)
library(reshape2)
library(scales)
library(dplyr)
library(gridExtra)
library(grid)
produce_all_plots = FALSE
```

In this section we will look into resistance patterns *Staphylococcus aureus*.


# Aggregating Data For Public Samples
First we'll get all publicly available *S. aureus* samples.

```{r}
ps <- get_public_samples()
```

# MRSA and MSSA
We defined MRSA by the presence of the *mecA*. Samples which did not have evidence for *mecA* were classified as MSSA.

## Primer based classification
First we'll use the results from the primer based SCCmec classification to identify samples with full matches to the *mecA* primers. It is important to note these results will only identify SCCmec types containing *mecA* (Example: SCCmec Xi has *mecC* and will not be included in these)

### Strict (Full Hits Only)
```{r}
sccmec_primer <- get_sccmec_type(ps$sample_id)
table(sccmec_primer$meca)
```

```{r}
sccmec_counts <- as.data.frame(colSums(sccmec_primer[,2:11]))
colnames(sccmec_counts) <- c('Total')
sccmec_counts <- data.frame(Type=rownames(sccmec_counts),
                            Total=sccmec_counts$Total)
sccmec_counts
```

### Relaxed (Hamming Distance)
```{r}
sccmec_type_hd <- get_sccmec_type(ps$sample_id, 
                                  hamming = TRUE)
table(sccmec_type_hd$meca)
```

## Protein Based Classification

```{r}
sccmec_proteins <- get_sccmec_protein_hits(ps$sample_id)
max_score <- group_by(sccmec_proteins,target) %>%
             summarise(maxscore = max(bitscore))
sccmec_proteins <- merge(sccmec_proteins, max_score,
                         by='target')
```

```{r}
sccmec_proteins$BSR <- sccmec_proteins$bitscore / sccmec_proteins$maxscore
table(sccmec_proteins[sccmec_proteins$BSR >0.95,]$target)
```

```{r fig.width=12, fig.asp=0.4}
p <- ggplot(data=sccmec_proteins[sccmec_proteins$BSR > 0.95],
            aes(x=target)) +
    ylab("Count") +
    xlab("SCCmec Target") +
    geom_bar() +
    theme_bw() +
    theme(axis.text=element_text(size=12),
          axis.title=element_text(size=14,face="bold"),
          legend.title = element_blank())
p
```

```{r fig.width=12, fig.asp=0.4}
a <- sccmec_proteins[sccmec_proteins$target =='mecA',]
sccmec_proteins_mec <- data.frame(sample_id=a$sample_id, BSR=a$BSR)
sccmec_proteins_mec$mec <- ifelse(
    sccmec_proteins_mec$BSR >= 0.95, TRUE, FALSE
)
p <- ggplot(data=sccmec_proteins_mec, aes(x=BSR)) +
    xlab("Blast Score Ratio") +
    ylab("Count") +
    geom_histogram(binwidth = 0.025) +
    theme_bw() +
    theme(axis.text=element_text(size=12),
          axis.title=element_text(size=14,face="bold"),
          legend.title = element_blank())
p
```

```{r fig.width=12, fig.asp=0.4}
a <- hist(sccmec_proteins_mec$BSR, plot=FALSE)
meca_bsr <- data.frame(
    region=sapply(1:length(a$counts), function(x){
        paste0(a$breaks[x], " - ",a$breaks[x+1])
    }), 
    count=a$counts
)
meca_bsr
```

## Ariba based classification
We can use the results from Ariba (via MEGARes) to identify samples with predicted resistance to methicillin. We will do this in two ways, first by only looking at results with a match (strict), and the other being those results that include partial assemblies (relaxed). A reminder, the Ariba results only include samples with paired end reads.

### Strict
These results are based on the a match to a SCCmec related cluster.
```{r}
ariba <- get_sccmec_ariba(ps$sample_id, resistance_report = TRUE)
table(ariba$mec)
```
### Relaxed
These results allow for partial matches to a SCCmec related cluster.
```{r}
ariba_relaxed <- get_sccmec_ariba(ps$sample_id,
                                  resistance_report = TRUE,
                                  include_all=TRUE)
table(ariba_relaxed$mec)
```

## SCCmec Cassette Coverage
```{r}
sccmec_coverage <- get_sccmec_cassette_coverages(ps$sample_id)
```

### Group By Most Covered SCCmec Type
```{r Table Top SCCmec Type Mapping}
top_type <- sccmec_coverage %>% group_by(sample_id) %>% slice(
    which.max(total)
)
table(top_type[top_type$total > 0.5,]$cassette)
```

```{r}
length(top_type[top_type$total > 0.5,]$cassette)
```

### Group By Most Covered *mec* Region
```{r}
top_mec <- sccmec_coverage %>% group_by(sample_id) %>% slice(
    which.max(meca_total)
)
table(top_mec[top_type$total > 0.5,]$cassette)
```

#### Plot Of Top SCCmec Covered and *mec* Region Covered

##### *mec* Predicted By Primers
```{r Figure SCCmec Primer, fig.width=12, fig.asp=0.4}
p <- ggplot(data=merge(top_type, sccmec_primer, by='sample_id'),
            aes(total, meca_total, colour = meca)) + 
        ylab("mec Region Covered") +
        xlab("Top SCCmec Type Covered") +
        geom_point() +
        scale_fill_manual(values=c("#2ca25f", "#5ab4ac")) +
        theme_bw() +
        theme(axis.text=element_text(size=12),
              axis.title=element_text(size=14,face="bold"),
              legend.title = element_blank())
p
```

##### *mec* Predicted By Ariba (Strict)
```{r Figure SCCmec Ariba, fig.width=12, fig.asp=0.4}
p <- ggplot(data=merge(top_type, ariba, by='sample_id'),
            aes(total, meca_total, colour = mec)) + 
        ylab("mec Region Covered") +
        xlab("Top SCCmec Type Covered") +
        geom_point() +
        scale_fill_manual(values=c("#2ca25f", "#5ab4ac")) +
        theme_bw() +
        theme(axis.text=element_text(size=12),
              axis.title=element_text(size=14,face="bold"),
              legend.title = element_blank())
p
```

##### *mec* Predicted By Ariba (Relaxed)
```{r fig.width=12, fig.asp=0.4}
p <- ggplot(data=merge(top_type, ariba_relaxed, by='sample_id'),
            aes(total, meca_total, colour = mec)) + 
        ylab("mec Region Covered") +
        xlab("Top SCCmec Type Covered") +
        geom_point() +
        scale_fill_manual(values=c("#2ca25f", "#5ab4ac")) +
        theme_bw() +
        theme(axis.text=element_text(size=12),
              axis.title=element_text(size=14,face="bold"),
              legend.title = element_blank())
p
```

#### Plot Of SCCmec Covered and Top *mec* Region Covered

##### *mec* Predicted By Primers
```{r fig.width=12, fig.asp=0.4}
p <- ggplot(data=merge(top_mec, sccmec_primer, by='sample_id'),
            aes(total, meca_total, colour = meca)) + 
        ylab("Top mec Region Covered") +
        xlab("SCCmec Type Covered") +
        geom_point() +
        scale_fill_manual(values=c("#2ca25f", "#5ab4ac")) +
        theme_bw() +
        theme(axis.text=element_text(size=12),
              axis.title=element_text(size=14,face="bold"),
              legend.title = element_blank())
p
```

##### *mec* Predicted By Ariba (Strict)
```{r fig.width=12, fig.asp=0.4}
p <- ggplot(data=merge(top_mec, ariba, by='sample_id'),
            aes(total, meca_total, colour = mec)) + 
        ylab("Top mec Region Covered") +
        xlab("SCCmec Type Covered") +
        geom_point() +
        scale_fill_manual(values=c("#2ca25f", "#5ab4ac")) +
        theme_bw() +
        theme(axis.text=element_text(size=12),
              axis.title=element_text(size=14,face="bold"),
              legend.title = element_blank())
p
```

##### *mec* Predicted By Ariba (Relaxed)
```{r fig.width=12, fig.asp=0.4}
p <- ggplot(data=merge(top_mec, ariba_relaxed, by='sample_id'),
            aes(total, meca_total, colour = mec)) + 
        ylab("Top mec Region Covered") +
        xlab("SCCmec Type Covered") +
        geom_point() +
        scale_fill_manual(values=c("#2ca25f", "#5ab4ac")) +
        theme_bw() +
        theme(axis.text=element_text(size=12),
              axis.title=element_text(size=14,face="bold"),
              legend.title = element_blank())
p
```

#### Plot of Each SCCmec Type Individually
##### Function For Plotting
```{r}
plot_by_sccmectype <- function(coverage, column) {
    p <- ggplot(data=coverage,
                aes(x=total, y=meca_total, colour = mec)) +
        ylab("mec Region Covered") +
        xlab(paste0("SCCmec Type ", column," Covered")) +
        geom_point() +
        scale_fill_manual(values=c("#2ca25f", "#5ab4ac")) +
        theme_bw() +
        theme(axis.text=element_text(size=12),
              axis.title=element_text(size=14,face="bold"),
              legend.title = element_blank())
    return(p)
}
```

#### Each SCCmec Type Individually (Colored By Ariba (Strict))
```{r fig.width=12, fig.asp=0.4}
if (produce_all_plots) {
    for (column in unique(sccmec_coverage$cassette)) {
        print(plot_by_sccmectype(
            merge(sccmec_coverage[sccmec_coverage$cassette == column,],
                  ariba, 
                  by='sample_id'),
            column
        ))
    }
}
```

#### Each SCCmec Type Individually (Colored By Ariba (Relaxed))
```{r fig.width=12, fig.asp=0.4}
if (produce_all_plots) {
    for (column in unique(sccmec_coverage$cassette)) {
        print(plot_by_sccmectype(
            merge(sccmec_coverage[sccmec_coverage$cassette == column,],
                  ariba_relaxed, 
                  by='sample_id'),
            column
        ))
    }
}
```

#### Compare *mec* Predictions
```{r}
mec <- merge(
    ps, 
    data.frame(
        sample_id=sccmec_proteins_mec$sample_id,
        protein_mec=sccmec_proteins_mec$mec
    ), 
    by='sample_id', 
    all=TRUE
)
mec[is.na(mec$protein_mec),]$protein_mec <- FALSE
mec <- merge(mec, data.frame(
    sample_id=sccmec_primer$sample_id, 
    primer_mec=sccmec_primer$meca
), by='sample_id', all=TRUE)
mec <- merge(mec, data.frame(
    sample_id=ariba$sample_id, 
    ariba_mec=ariba$mec
), by='sample_id', all=TRUE)
mec$agreement <- paste0(
    ifelse(mec$protein_mec, 1, 0),
    ifelse(mec$primer_mec, 1, 0),
    ifelse(is.na(mec$ariba_mec), '-', ifelse(
        mec$ariba_mec, 1, 0
    ))
)
table(mec$agreement)
```

*Notes*

* 0: **mec** not predicted
* 1: **mec** predicted
* -: Not tested by Ariba (Single-End reads)

The order of numbering is:

* 1: **mecA** with BSR > 0.95 based on Proteins
* 2: **mecA** with full Primer hit
* 3: **mec** with full match based on Ariba

Example:

* 00-- : Single-End, Protein and Primer are False
* 0000 : All approaches agree that **mec** is not predicted
* 1111 : All approaches agree that **mec** is predicted

```{r}
metrics <- merge(
    ps,
    merge(
        get_assembly_stats(ps$sample_id),
        merge(
            get_metadata(ps$sample_id),
            get_sequence_quality(ps$sample_id, stage='cleanup'),
            by='sample_id'
        ),
        by='sample_id'
    ),
    by='sample_id'
)
metrics$year <- sapply(
    metrics$first_public,
    function(x) {
        strsplit(x, "-")[[1]][1]
    }
)

metrics$rank_name <- ifelse(
    metrics$rank.x == 3,
    'Gold',
    ifelse(
        metrics$rank.x == 2,
        'Silver',
        'Bronze'
    )
)
platform <- metrics[,c('sample_id', 'instrument_model', 'total_contig', 'study_accession', 'rank_name')]
mec_temp <- merge(mec, platform, by='sample_id')
```

```{r}
table(mec_temp[mec_temp$agreement == '111',]$study_accession)
```

```{r}
summary(mec_temp[mec_temp$agreement == '001',]$total_contig)
```

```{r}
colnames(metrics)
```

```{r}
mec <- merge(mec, data.frame(
    sample_id=top_type$sample_id, 
    total=top_type$total
), by='sample_id')
mec <- merge(mec, data.frame(
    sample_id=top_type$sample_id, 
    mec_total=top_type$meca_total
), by='sample_id')
table(mec[mec$total >= 0.5,]$primer_mec)
```
```{r}
table(mec[mec$total >= 0.5,]$protein_mec)
```

```{r}
table(mec[mec$total >= 0.5,]$ariba_mec)
```

```{r}
table(mec[mec$total >= 0.5,]$ariba_relaxed_mec)
```

SCCmec cassettes with 50% coverage but 0% in *mec* region.
```{r}
table(mec[mec$total >= 0.5 & mec$mec_total == 0,]$ariba_mec)
```

#### *mecA* Presence By ST
```{r}
meca_groups <- merge(
    data.frame(
        sample_id=ps$sample_id, 
        st=ps$st, 
        rank=ps$rank
    ),
    ariba,
    by='sample_id'
)
meca_groups$status <- ifelse(
    meca_groups$mec == TRUE, 'MRSA', 'MSSA'
)
meca_by_st <- plyr::count(meca_groups, c('st', 'status'))
head(meca_by_st)
```

#### Get the Top 10 STs
```{r}
top_st <- get_top_sequence_types()
top_st
```

#### MRSA/MSSA Distribution For Top 10 Sequence Types
Now we are ready to plot out the distribution of MRSA and MSSA predictions by sequence type.

```{r fig.width=12, fig.asp=0.4}
top_st_meca <- merge(top_st, meca_by_st, by='st')
p <- ggplot(data=top_st_meca, aes(x=reorder(st, -count), y=freq,
                                  fill = status)) +
    xlab("Sequence Type") +
    ylab("Count") +
    geom_bar(stat="identity", position = "dodge") +
    geom_text(aes(label=freq), vjust = -0.5,
              position = position_dodge(.9)) + 
    scale_fill_manual(values=c("#2ca25f", "#5ab4ac")) +
    theme_bw() +
    theme(axis.text=element_text(size=12),
          axis.title=element_text(size=14,face="bold"),
          legend.title = element_blank())
p
```

#### MRSA/MSSA Distribution For Top 10 Sequence Types
Now we are ready to plot out the distribution of MRSA and MSSA predictions by sequence type.

```{r fig.width=8, fig.asp=0.4}
top_st_meca <- merge(top_st, meca_by_st, by='st')
counts <- top_st_meca[top_st_meca$status == 'MRSA',]
counts$label <- paste0(counts$freq, ", ", counts$count - counts$freq)
counts <- counts[,c('st', 'count', 'label')]
p <- ggplot(data=top_st_meca, aes(x=reorder(st, -count), y=freq,
                                  fill = status, label=freq)) +
    xlab("Sequence Type") +
    ylab("Count") +
    geom_bar(stat="identity", position = "stack") +
    geom_text(size=3, aes(as.character(st), count + 200, label = label, fill = NULL), data = counts) +
    scale_fill_manual(values=c("#2ca25f", "#5ab4ac")) +
    theme_bw() +
    theme(axis.text=element_text(size=12),
          axis.title=element_text(size=14,face="bold"),
          legend.title = element_blank())
p
```

### Resistance Patterns Based on Ariba and MegaRes

#### Resistance Classes With A Match
```{r}
resistance_report <- get_resistance_results(ps$sample_id,
                                            resistance_report=TRUE)
df <- as.data.frame(colSums(
    resistance_report[,2:ncol(resistance_report)]
))
colnames(df) <- c('total')
resistance_report_counts <- data.frame(
    class=rownames(df),
    total=df$total
)
resistance_report_counts[resistance_report_counts$total > 0,]
```

#### Resistance Classes Without A Match
```{r}
data.frame(class=as.character(
    resistance_report_counts[resistance_report_counts$total == 0,]$class
))
```

```{r}
cluster_report <- get_resistance_results(ps$sample_id,
                                         cluster_report=TRUE)
df <- as.data.frame(colSums(
    cluster_report[,2:ncol(cluster_report)]
))
colnames(df) <- c('total')
cluster_report_counts <- data.frame(
    cluster=rownames(df),
    total=df$total
)
cluster_report_counts[cluster_report_counts$total > 0,]
```

#### Resistance Clusters Without A Match
```{r}
data.frame(class=as.character(
    cluster_report_counts[cluster_report_counts$total == 0,]$cluster
))
```


#### Group
```{r}
resistance_groups <- merge(
    data.frame(sample_id=ps$sample_id, st=ps$st, rank=ps$rank),
    resistance_report,
    by='sample_id'
)
```

#### By ST

#### Function For Plotting
```{r}
plot_by_st <- function(group, top_st, column) {
    by_st <- plyr::count(
        group, 
        c('st',
          ifelse(
              length(strsplit(column, ' ')[[1]]) >= 2, 
              paste0("`", column, "`"), 
              column
          )
        )
    )
    by_st$status <- ifelse(
        by_st[,make.names(column)] == TRUE, 'Resistant', 'Susceptible'
    )
    top_st_resistance <- merge(top_st, by_st, by='st')
    p <- ggplot(data=top_st_resistance, aes(x=reorder(st, -count),
                                            y=freq, fill = status)) +
        xlab("Sequence Type") +
        ylab(paste0("Count")) +
        geom_bar(stat="identity", position = "stack") +
        scale_fill_manual(values=c("#2ca25f", "#5ab4ac"),
                          name = column) +
        theme_bw() +
        theme(axis.text=element_text(size=12),
              axis.title=element_text(size=14, face="bold"),
              legend.position="top",
              legend.title=element_blank())
    return(p)
}
```

#### Resistance Classes With A Match By Top 10 Sequence Types
```{r fig.width=12, fig.asp=0.4}
top_st <- get_top_sequence_types()
for (column in sort(as.character(
        resistance_report_counts[resistance_report_counts$total > 0,]$class
    ))) {
    p <- plot_by_st(resistance_groups, top_st, column)
    print(column)
    print(p)
}
```

```{r figure 5, fig.width=12, fig.asp=0.8}
grid_arrange_shared <- function(..., ncol = length(list(...)), nrow = 1, position = c("bottom", "right")) {
  plots <- list(...)
  position <- match.arg(position)
  g <- ggplotGrob(plots[[1]] + theme(legend.position = position))$grobs
  legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
  lheight <- sum(legend$height)
  lwidth <- sum(legend$width)
  gl <- lapply(plots, function(x) x + theme(legend.position="none"))
  gl[[1]] <-arrangeGrob(gl[[1]], top = textGrob(
        "A", x = unit(0, "npc"), y = unit(1, "npc"), just=c("left","top"),
         gp=gpar(col="black", fontsize=18, fontface="bold"))
    )
  gl[[2]] <-arrangeGrob(gl[[2]], top = textGrob(
        "B", x = unit(0, "npc"), y = unit(1, "npc"), just=c("left","top"),
         gp=gpar(col="black", fontsize=18, fontface="bold"))
    )
  gl[[3]] <-arrangeGrob(gl[[3]], top = textGrob(
        "C", x = unit(0, "npc"), y = unit(1, "npc"), just=c("left","top"),
         gp=gpar(col="black", fontsize=18, fontface="bold"))
    )
  gl[[4]] <-arrangeGrob(gl[[4]], top = textGrob(
        "D", x = unit(0, "npc"), y = unit(1, "npc"), just=c("left","top"),
         gp=gpar(col="black", fontsize=18, fontface="bold"))
    )
      
  gl <- c(gl, ncol = ncol, nrow = nrow)

  combined <- switch(position,
                     "bottom" = arrangeGrob(do.call(arrangeGrob, gl),
                                            legend,
                                            ncol = 1,
                                            heights = unit.c(unit(1, "npc") - lheight, lheight)),
                     "right" = arrangeGrob(do.call(arrangeGrob, gl),
                                           legend,
                                           ncol = 2,
                                           widths = unit.c(unit(1, "npc") - lwidth, lwidth)))

  grid.newpage()
  grid.draw(combined)

  # return gtable invisibly
  invisible(combined)
}

top_st <- get_top_sequence_types()
p1 <- plot_by_st(resistance_groups, top_st, 'Betalactams (mec)')
p2 <- plot_by_st(resistance_groups, top_st, 'Aminoglycosides')
p3 <- plot_by_st(resistance_groups, top_st, 'Fosfomycin')
p4 <- plot_by_st(resistance_groups, top_st, 'MLS')
grid_arrange_shared(p1,p2,p3,p4, nrow=2, ncol=2)
# Output plot to PDF and PNG
pdf(paste0(
    getwd(),
    '/../figures/figure-04-resistance-facet-top-10-sequence-types.pdf'
), height=6, width=12)
grid_arrange_shared(p1,p2,p3,p4, nrow=2, ncol=2)
dev_null <- dev.off()

png(paste0(
    getwd(),
    '/../figures/figure-04-resistance-facet-top-10-sequence-types.png'
), height=600, width=1200)
grid_arrange_shared(p1,p2,p3,p4, nrow=2, ncol=2)
dev_null <- dev.off()
```

## Session Info
```{r}
sessionInfo()
```
---
title: "Results Section: Public Metadata"
output:
  pdf_document: default
---
```{r}
library(staphopia)
library(ggplot2)
library(reshape2)
```

## Aggregating Data For Public Samples
First we'll get all publicly available *S. aureus* samples.

```{r}
ps <- get_public_samples()
```

## Variation From *S. aureus* N315
In Staphopia all samples had variants (SNPs and InDels) called using *S. aureus* N315 as the reference genome. In this section we'll visualize the total number of variants each sample has. This will give us an idea of the sequenced genitic diversity with respect to N315.

### Gather Variant Counts
We will use `get_variant_counts()` to get the variant counts for each sample. We will also order the counts by the total.

```{r}
variant_counts <- get_variant_counts(ps$sample_id)
variant_counts <- variant_counts[order(total),]
```

### Summary of Variant Counts
#### Total Variants (SNPs and InDels)
```{r}
summary(variant_counts$total)
```

#### SNPs
```{r}
summary(variant_counts$snp_count)
```

#### InDels
```{r}
summary(variant_counts$indel_count)
```

### Visualizing Variant Counts
#### Total Variants (SNPs and InDels)
```{r fig.width=12, fig.asp=0.4}
p <- ggplot(data=variant_counts, aes(x=seq(1,nrow(variant_counts)), y=total)) +
    xlab("Sample Count") +
    ylab("Count") +
    geom_bar(stat='identity') +
    scale_x_continuous(breaks = seq(0, nrow(variant_counts), by = 5000)) +
    scale_y_continuous(breaks = seq(0, max(variant_counts$total), by=10000),
                       labels = scales::comma) +
    theme_bw() +
    theme(axis.text=element_text(size=12),
          axis.title=element_text(size=14,face="bold"))
p
```

#### SNPs Only
```{r figure 9, fig.width=12, fig.asp=0.4}
cutoff <- max(variant_counts[variant_counts$snp_count < 60000,]$snp_count)
variant_counts$fill <- ifelse(variant_counts$snp_count > cutoff, TRUE, FALSE)
p <- ggplot(data=variant_counts, aes(x=seq(1,nrow(variant_counts)), y=snp_count,
                                     fill=fill)) +
    xlab("Samples") +
    ylab("Count") +
    geom_hline(yintercept = cutoff, linetype="dotted") +
    geom_bar(stat='identity') +
    annotate("text", x = 40000, y = 60000, label = "S. argenteus", fontface=3) +
    annotate("text", x = 20000, y = 50000, label = "S. aureus", fontface=3) +
    scale_x_continuous(breaks = seq(0, nrow(variant_counts), by = 5000)) +
    scale_y_continuous(breaks = seq(0, max(variant_counts$snp_count), by=20000),
                       labels = scales::comma) +
    scale_fill_manual(values=c("#9C8443", "#B9C6C6")) +
    theme_bw() +
    theme(axis.text=element_text(size=12),
          axis.title=element_text(size=14,face="bold"),
          legend.position="none")
p

# Output plot to PDF and PNG
staphopia::write_plot(p, paste0(getwd(), '/../figures/figure-05-snp-accumulation'))
```

```{r fig.width=12, fig.asp=0.4}
cutoff <- max(variant_counts[variant_counts$snp_count < 60000,]$snp_count)
variant_counts$fill <- ifelse(variant_counts$snp_count > cutoff, TRUE, FALSE)
p <- ggplot(data=variant_counts[variant_counts$snp_count > 52000,], aes(
        x=seq(1,nrow(variant_counts[variant_counts$snp_count > 52000,])),
        y=snp_count,
        fill=fill)
    ) +
    xlab("Samples") +
    ylab("Count") +
    geom_hline(yintercept = cutoff, linetype="dotted") +
    geom_bar(stat='identity') +
    scale_x_continuous(breaks = seq(
        0, nrow(variant_counts[variant_counts$snp_count > 52000,]), by = 20)
    ) +
    scale_y_continuous(breaks = seq(0, max(variant_counts$snp_count), by=10000),
                       labels = scales::comma) +
    scale_fill_manual(values=c("#9C8443", "#B9C6C6")) +
    theme_bw() +
    theme(axis.text=element_text(size=12),
          axis.title=element_text(size=14,face="bold"),
          legend.position="none")
p
```

#### InDels Only
```{r fig.width=12, fig.asp=0.4}
p <- ggplot(data=variant_counts, aes(x=seq(1,nrow(variant_counts)), y=indel_count)) +
    xlab("Sample Count") +
    ylab("Count") +
    geom_bar(stat='identity') +
    scale_x_continuous(breaks = seq(0, nrow(variant_counts), by = 5000)) +
    scale_y_continuous(breaks = seq(0, max(variant_counts$indel_count), by=500),
                       labels = scales::comma) +
    theme_bw() +
    theme(axis.text=element_text(size=12),
          axis.title=element_text(size=14,face="bold"))
p
```

# Non-Redundant Diversity Dataset
Below are the set of commands used to generate the NRD dataset. These were executed on a local server using Staphopia functions not available through the API due to execution times.

```{bash NRD Dataset, eval=FALSE}
# Create NRD dataset (non-redundant-set.txt)
~/staphopia-web/manage.py create_nonredundant_set nrd-dataset --published --rank 3 \
                                                              --singletons

# Extract all genes without InDels and validate SNPs via 31-mers
NRD="nrd-annotation-summary.txt"
STAPHOPIA="/home/rpetit/cgc-staphopia/staphopia"
mkdir extract-genome
awk '{print $8}' non-redundant-set.txt | \
    grep -v sample_id | \
    xargs -I {} ~/staphopia-web/manage.py extract_genome {} $NRD $STAPHOPIA

# Identify genes with less than 3 samples with missing 31-mers
cat extract-genome/*-annotation-details.txt | \
    awk '{if ($3=="False"){print $2}}' | \
    sort | \
    uniq -c | \
    awk '{if ($1 >= 377){print $2}}' | \
    sort -n > validated-genes.txt

# Extract validated genes (extract-core.fasta, extract-core.phylip)
~/staphopia-web/manage.py extract_core validated-genes.txt extract-genome/

# Build Phylogeny
# Start Tree
iqtree -s extract-core.phylip -m GTR -nt AUTO -fast -pre start-tree | \
    tee iqtree-start.stdout.txt

# Identify Recombination
ClonalFrameML start-tree.treefile extract-core.fasta clonalframe -emsim 100 | \
    tee clonalframe.stdout.txt

# Remove Recombination
maskrc-svg.py clonalframe --aln extract-core.fasta \
                          --out extract-core-cfmasked.fasta --symbol '-'| \
    tee maskrc-svg.stdout.txt

# Final Tree
iqtree -s extract-core-cfmasked.fasta -m GTR -nt 11 -pre final-tree \
       -bb 1000 -alrt 1000 -wbt -wbtl -alninfo | \
    tee iqtree-final.stdout.txt
```

# Session Info
```{r}
sessionInfo()
```
---
title: "Results Section: Pipeline Design and Processing 43,000+ genomes"
output:
  pdf_document: default
---

In this notebook will be generating statistics and plots related to processing 43,000+ 
genomes on Seven Bridges Cancer Genomics Cloud (CGC) platform.

### Load Up Packages
```{r}
library(staphopia)
library(ggplot2)
library(dplyr)
```

### Read In The Data
```{r}
results <- read.table("../data/cgc-runs.txt", header = TRUE, sep = "\t")
colnames(results)
```

This leaves use with 9 columns:

1. name: Name of the job
2. status: Job's status
3. project: CGC project job was executed from.
4. app: CGC app used to execute the job.
5. created_by: User who submitted the job.
6. total_time: Total amount of time (in minutes) a job was queued and run
7. run_time: Total amount of time (in minutes) a job took to complete
8. queue_time: Total amount of time (in minutes) a job was queued
9. price: Total cost of the run

### Clean Up The Data
Before we generate statistics and plots, we need to clean the data. There are jobs where the *run_time* 
and *price* were not properly reported from CGC. We will filter samples where the *run_time* is 0.

```{r}
results_clean <- results[results$run_time > 0, ]
nrow(results) - nrow(results_clean)
```

### Job Summary
#### Run Time Summary
```{r}
summary(results_clean$run_time)
```

#### Number of Jobs With > 120 Minute Runtime
```{r}
nrow(results_clean[results_clean$run_time > 120, ])
```

#### Summary of Jobs With Run Time Between 10 and 120 Minutes
```{r}
summary(results_clean[between(results_clean$run_time, 10, 120), ]$run_time)
```

### Plots
#### Run Time (Complete)
```{r fig.width=12, fig.asp=0.4}
p <- ggplot(data=results_clean, aes(run_time)) +
    xlab("Run Time (In Minutes)") +
    ylab("Count") +
    geom_histogram(bins=100) +
    theme_bw() +
    theme(axis.text=element_text(size=12),
          axis.title=element_text(size=14,face="bold"))
p
```

#### Pipeline Run Time (Between 10-120 Minutes)
```{r Figure 2 Pipeline Run Time, fig.width=12, fig.asp=0.4}
p <- ggplot(data=results_clean[between(results_clean$run_time, 10, 120),],
            aes(run_time)) +
    xlab("Run Time (In Minutes)") +
    ylab("Count") +
    geom_histogram(bins=100) +
    theme_bw() +
    theme(axis.text=element_text(size=12),
          axis.title=element_text(size=14,face="bold"))
p

# Output plot to PDF and PNG
staphopia::write_plot(p, paste0(getwd(), '/../figures/supplementary-figure-02-pipeline-run-time'))
```


## Session Info
```{r}
sessionInfo()
```
