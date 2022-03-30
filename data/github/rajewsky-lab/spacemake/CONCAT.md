<a href="https://pypi.org/project/spacemake/">
   <img src="https://img.shields.io/pypi/v/spacemake.svg" / ></a>
   
<a href="https://spacemake.readthedocs.io/">
   <img src="https://readthedocs.org/projects/spacemake/badge/?version=latest" / ></a>
   
 <a href="https://pepy.tech/project/spacemake">
   <img src="https://pepy.tech/badge/spacemake" / ></a>

# spacemake: processing and analysis of large-scale spatial transcriptomics data

<img src="https://raw.githubusercontent.com/rajewsky-lab/spacemake/master/docs/graphical_abstract_twitter.png" width="500" />

Documentation can be found [here](https://spacemake.readthedocs.io/en/latest/).

---
output:
  html_document:
    toc: true
    self_contained: yes
version: 0.1.1
author: Tamas Ryszard Sztanka-Toth, Nikolaos Karaiskos
email: tamasryszard.sztanka-toth@mdc-berlin.de, nikolaos.karaiskos@mdc.berlin.de
license: GPL
---

```{r knitr_options, include=FALSE, cache=FALSE}
knitr::opts_chunk$set(
  autodep = TRUE,
  message = FALSE,
  warning = FALSE,
  comment = NA
)
```

```{r libraries, include = F, cache=F}
library(tidyverse)
library(yaml)
library(knitr)
library(magrittr)
library(kableExtra)
library(cowplot)

theme_set(theme_cowplot(18))

cpalette <- list('orange' = '#D55E00', 'blue' = '#0072B2', 'green' = '#009E73', 'black' = '#000000', 'yellow' = '#F0E442', 
				 'grey' = '#999999', 'light_orange' = "#E69F00", 'light_blue' = "#56B4E9")

clrs = c('umis'=cpalette$light_orange, 'genes' = cpalette$light_blue,
         'reads'=cpalette$green, 'pcr'=cpalette$pink, 'pct_counts_mt'= 'black')

source(snakemake@params$r_shared_scripts)
```

```{r read_data, echo = F}
obs_df <- read_csv(snakemake@input$obs_df) %>%
    mutate(pct_counts_mt = as.double(pct_counts_mt))
var_df <- read_csv(snakemake@input$var_df)

n_beads = nrow(obs_df)
data_empty <- n_beads == 0

if (!data_empty){
    data_empty <- sum(obs_df$total_counts) < 1
}
```

### Overview

```{r create_overview_table, echo = F}
puck_width_um <- snakemake@params$puck_variables$width_um
parameter_stats <- rbind(c('UMI filter', snakemake@wildcards$umi_cutoff),
          c('# genes in data', nrow(var_df)),
          c('# of spots in data', obs_df %>% select(cell_bc) %>% unique %>% nrow),
          c('median UMI', median(obs_df$total_counts)),
          c('median genes', median(obs_df$n_genes_by_counts)),
          c('puck width (um)', puck_width_um))

parameter_stats %>%
    kbl(col.names=NULL) %>%
    kable_classic_2(full_width=F) %>%
    #pack_rows('Sample info', 1, 7) %>%
    pack_rows('Data info', 1, 5)
```

### Histogram of metrics over beads

```{r data_empty, echo = F, eval = data_empty, results='asis'}
cat('This dataset has 0 beads passing the filters')
```

```{r plot_histogram_of_metrics, echo = F, fig.width=10, eval=!data_empty}
pl1 <- obs_df %>%
    select(cell_bc, total_counts, n_genes_by_counts) %>%
    dplyr::rename(umis=total_counts,
                  genes=n_genes_by_counts) %>%
    gather('obs', 'value', -cell_bc) %>%
    ggplot(aes(value, fill=obs)) +
        geom_histogram(bins=100) +
        scale_fill_manual(values=clrs) +
        facet_wrap(~obs, ncol=2, scales='free')

if (length(unique(obs_df$total_counts)) > 50){
    pl1 <- pl1 +
        scale_x_log10() +
        annotation_logticks(sides='b') +
        labs(x='log(value)', y='count')
} else {
    pl1 <- pl1 +
        labs(x ='value', y='count')
}

pl2 <- obs_df %>%
    select(cell_bc, pct_counts_mt) %>%
    gather('obs', 'value', -cell_bc) %>%
    ggplot(aes(value, fill=obs)) +
        geom_histogram(bins=100) +
        scale_fill_manual(values=clrs) +
        facet_wrap(~obs, ncol=2, scales='free')+
        labs(x='value', y='count')

plot_grid(pl1,
          plot_grid(pl2, NULL, ncol=2, rel_widths=c(1.5,1)),
          ncol=1, rel_heights=c(1,1), labels="")
```


```{r check_completeness, echo = F, fig.size=10, fig.width=10}
obs_df <- obs_df %>%
    gather('res', 'cluster', starts_with('leiden'))

data_complete = F

if ('cluster' %in% colnames(obs_df)){
    data_complete = T
}

is_spatial <- snakemake@params$is_spatial & data_complete

fig_height=6
fig_width=7
```

```{r incompleteness_warning, eval=!data_complete, echo=F, results='asis'}
cat('### WARNING: data incomplete\n\n')
cat(paste0('This dataset has ', n_beads, ' cells passing the filter of ', snakemake@wildcards$umi_cutoff, ' and ', nrow(var_df), ' genes.\n\n'))
cat(paste0('This dataset is too small, so it couldn\'t be properly clustered and analysed automatically'))
```


```{r umap_text, eval=data_complete, echo =F, results='asis'}
cat('### UMAP plots using different resolutions {.tabset}\n\n')
```


```{r plot_umi_umap, eval=data_complete, echo=F, fig.height=fig_height, fig.width=fig_width}
obs_df <- obs_df %>%
    mutate(cluster = factor(cluster)) %>%
    mutate(log1p_total_counts = log2(1+total_counts))

n_cells <- obs_df %>%
    select(cell_bc) %>%
    unique %>%
    nrow()

def_plot_bead_size <- ifelse(n_cells > 5000, 0.4, 0.75)
def_plot_bead_size <- ifelse(n_cells > 10000, 0.1, def_plot_bead_size)
def_plot_bead_size <- ifelse(n_cells > 25000, 0.05, def_plot_bead_size)

obs_colnames <- colnames(obs_df)

# barcode file attached at the python level
is_hexagonal <- snakemake@params$run_mode_variables$mesh_type == 'hexagon' & snakemake@params$run_mode_variables$mesh_data

filter_quant <- quantile(obs_df$total_counts, 0.9, na.rm=T)
limits <- c(0, quantile(obs_df$total_counts, 0.91, na.rm=T))

if (is_spatial){

    # calculate breaks and limits for the puck
    puck_file <- read_csv(snakemake@input$puck_file)
    x_limits <- puck_file$x_pos %>% {c(min(.), max(.))}
    y_limits <- puck_file$y_pos %>% {c(min(.), max(.))}

    ratio <- (x_limits[2] - x_limits[1] ) / (y_limits[2] - y_limits[1])

    scale_factor <- ifelse(puck_width_um < 3000, 2, 3)
    mm_dist <- max(10^scale_factor, round(puck_width_um/3, digits =-scale_factor))
    mm_diff <- mm_dist / 1000

    x_mm_breaks <- seq(0, puck_width_um, mm_dist)
    x_mm_breaks <- paste0(x_mm_breaks * mm_diff / mm_dist, 'mm')
    y_mm_breaks <- seq(0, puck_width_um / ratio, mm_dist)
    y_mm_breaks <- paste0(y_mm_breaks * mm_diff / mm_dist, 'mm')

    px_to_um <- (x_limits[2] - x_limits[1]) / snakemake@params$puck_variables$width_um

    x_breaks <- seq(x_limits[1], x_limits[2], px_to_um * mm_dist)
    y_breaks <- seq(y_limits[1], y_limits[2], px_to_um * mm_dist)

    puck_bead_size <- min(def_plot_bead_size, ifelse(snakemake@params$run_mode_variables$mesh_data,
        snakemake@params$run_mode_variables$mesh_spot_diameter_um / 40,
        snakemake@params$puck_variables$spot_diameter_um / 40))
    
    res_colnames <- obs_df$res %>%
        unique

    umi_pl <- obs_df %>%
        filter(res == res_colnames[1]) %>%
        arrange(total_counts) %>%
        mutate(total_counts = ifelse(total_counts > filter_quant, filter_quant, total_counts))

    if(is_hexagonal){
        umi_pl <- umi_pl %>%
            ggplot(aes(x=x_pos, y=y_pos, fill = total_counts, group=1)) +
                geom_hex(stat='identity', color='gray95') +
                coord_fixed()+
                scale_x_continuous(labels = x_mm_breaks, breaks = x_breaks, limits=x_limits) +
                scale_y_continuous(labels = y_mm_breaks, breaks = y_breaks, limits=y_limits) +
                scale_fill_viridis_c(option =  "magma", limits = limits) +
                guides(fill = guide_colorbar(barheight = 15)) + 
                labs(fill='UMI count', x='', y='')

    }
    else{
        umi_pl <- umi_pl %>%
            ggplot(aes(x_pos, y_pos, color = total_counts)) +
                geom_point(size=puck_bead_size) + 
                coord_fixed()+
                scale_x_continuous(labels = x_mm_breaks, breaks = x_breaks, limits=x_limits) +
                scale_y_continuous(labels = y_mm_breaks, breaks = y_breaks, limits=y_limits) +
                scale_color_viridis_c(option =  "magma", limits = limits) +
                guides(color = guide_colorbar(barheight = 15)) + 
                labs(color='# of UMIs\nper spatial unit', x='', y='')
    }
    umi_pl <- umi_pl + 
        theme(panel.background = element_rect(fill = 'gray95'), 
             legend.spacing = unit(0.1, 'cm'),
            axis.line = element_line(color = 'black'),
            text = element_text(color='black', size=18))
    umi_pl
}

```

```{r plot_clusters_umap_puck, echo =F, fig.height=fig_height, fig.width=fig_width, eval=data_complete, results='asis'}
library(pals)

cluster_clrs <- unname(glasbey())

#top10_marker_table <- read_table2(snakemake@input$cluster_markers)

for (i in obs_df %$% res %>% unique){
    res <- as.double(strsplit(i, '_')[[1]][2])
    cat(paste0('\n\n#### ', res, ' resolution {.tabset}\n\n'))
    dat <- obs_df %>%
        filter(res == i) %>%
        dplyr::select(cell_bc, umap_0, umap_1, cluster)

    umap_plot <- dat %>%
        ggplot(aes(umap_0, umap_1, color = cluster)) +
            geom_point(size=def_plot_bead_size) +
            guides(colour = guide_legend(override.aes = list(size=3)))+
            coord_fixed() +
            theme(axis.text = element_blank(),
                  legend.position = 'none',
                axis.ticks = element_blank(), axis.line = element_line(color='white'))

    n_clusters <- length(unique(dat$cluster))
    if(n_clusters< length(cluster_clrs)){
        umap_plot <- umap_plot + scale_color_manual(values=cluster_clrs)
    }

    if (is_spatial) {
        physical_plot <- obs_df %>%
            filter(res == i) %>%
            dplyr::select(cell_bc, x_pos, y_pos, cluster)

        if(is_hexagonal){
            physical_plot <- physical_plot %>%
            ggplot(aes(x_pos, y_pos, fill = cluster)) +
                geom_hex(stat='identity', color='gray95') +
                guides(fill = guide_legend(override.aes = list(size=3), ncol=2))+
                coord_fixed()+
                labs(x='', y='') +
                scale_x_continuous(labels = x_mm_breaks, breaks = x_breaks, limits=x_limits) +
                scale_y_continuous(labels = y_mm_breaks, breaks = y_breaks, limits=y_limits)
            if(n_clusters< length(cluster_clrs)){
                physical_plot <- physical_plot + scale_fill_manual(values=cluster_clrs)
            }
        } else{
            physical_plot <- physical_plot %>%
            ggplot(aes(x_pos, y_pos, color = cluster)) +
                geom_point(size=puck_bead_size) +
                guides(colour = guide_legend(override.aes = list(size=3), ncol=2))+
                coord_fixed()+
                labs(x='', y='') +
                scale_x_continuous(labels = x_mm_breaks, breaks = x_breaks, limits=x_limits) +
                scale_y_continuous(labels = y_mm_breaks, breaks = y_breaks, limits=y_limits)
            if(n_clusters< length(cluster_clrs)){
                physical_plot <- physical_plot + scale_color_manual(values=cluster_clrs)
            }

        }
        physical_plot <- physical_plot +
            theme(panel.background = element_rect(fill = 'gray95'),
                legend.spacing = unit(0.1, 'cm'),
                axis.line = element_line(color = 'black'),
                plot.subtitle=element_text(size=18),
                text = element_text(color='black')) +
            ggtitle('', subtitle=paste0('resolution = ', res))

        print(umap_plot)
        print(physical_plot)
    } else{
        print(umap_plot)
    }
}
```

```{r nhood_enrichment_title, echo =F, eval=is_spatial, results='asis'}
cat('\n\n### Neighborhood enrichment\n\n')
```

```{r plot_nhood_enrichment, echo =F, eval=is_spatial}
nhood_dat <- read_csv(snakemake@input$nhood_enrichment)

for (i in obs_df %$% res %>% unique){
    res <- as.double(strsplit(i, '_')[[1]][2])
    dat <- nhood_dat %>%
        filter(resolution == res) %>%
        mutate(a = as.character(1000 + cluster_a),
               b = as.character(1000 + cluster_b)) %>%
        mutate(zscore = ifelse(zscore > 100, 100, zscore),
               zscore = ifelse(zscore < -50, -50, zscore))

    labs <- dat %>%
        select(cluster_a, cluster_b, a, b) %>%
        filter(cluster_a %% 5 == 0, cluster_b %% 5 == 0)

    pl <- dat %>%
        ggplot(aes(a, b, fill=zscore)) +
            geom_tile(color='white') + 
            scale_fill_viridis_c(option='inferno', limits = c(-51, 101)) +
        guides(fill = guide_colorbar(barheight = 12)) + 
        scale_x_discrete(labels = labs$cluster_a, breaks=labs$a) +
        scale_y_discrete(labels = labs$cluster_b, breaks=labs$b) +
        coord_fixed() +
        labs(fill='neighborhood\nenrichment\nscore', x='cluster identity', y='cluster identity') +
        ggtitle('', subtitle=paste0('resolution = ', res))+
        coord_equal()+
        theme(plot.subtitle=element_text(size=18))
    print(pl)
}
```

---
output:
    html_document:
        toc: true
        toc_depth: 6
version: 0.2
author: Tamas Ryszard Sztanka-Toth, Nikolaos Karaiskos
email: tamasryszard.sztanka-toth@mdc-berlin.de, nikolaos.karaiskos@mdc.berlin.de
license: GPL
---

```{r knitr_options, include=FALSE, cache=FALSE}
knitr::opts_chunk$set(
  autodep = TRUE,
  message = FALSE,
  warning = FALSE,
  comment = NA
)
```

```{r load_libraries, echo =F}
library(tidyverse)
library(yaml)
library(cowplot)
library(knitr)
library(grid)
library(gtable)
library(kableExtra)
library(pals)

theme_set(theme_cowplot(18))
```

```{r functions, echo = F}
readStarLog <- function(log_file){
		out = list()
		lines = readLines(log_file)
	
		out$input_reads = (lines[6] %>% strsplit('\t') %>% unlist)[2] %>% as.integer

		out$uniq_mapped_reads = (lines[9] %>% strsplit('\t') %>% unlist)[2] %>% as.integer

		#out$avg_length = (lines[11] %>% strsplit('\t') %>% unlist)[2] %>% as.numeric
		
        tibble(observation=names(out), value=unlist(unname(out)))
	}
```

```{r load_run_info, echo = F}
run_mode_settings <- bind_rows(snakemake@params$run_mode) %>%
    nest(umi_cutoff = c(umi_cutoff)) %>%
    bind_cols(tibble(run_mode = names(snakemake@params$run_mode))) %>%
    gather('variable_name', 'value', -run_mode) %>%
    spread(run_mode, value)

run_mode_names <- run_mode_settings[-1] %>% names

sample_info <- bind_rows(snakemake@params$sample_info) %>%
    select(species, sequencing_date, investigator, experiment) %>%
    unique %>%
    gather('info_name', 'value') %>%
    add_row(info_name = 'project_id',
            value = snakemake@wildcards$project, .before=1)%>%
    add_row(info_name = 'sample_id',
            value = snakemake@wildcards$sample, .after=1)%>%
    add_row(info_name = 'puck_id',
            value = snakemake@wildcards$puck, .after=2)
```


## Run information

```{r print_sample_info, echo = F}
sample_info %>%
    kbl(col.names=NULL) %>%
    kable_classic_2(full_width=F, position='left')
```

saturation analysis v.0.1.1, generated on `r format(Sys.time(), '%d/%B/%Y')`

contact: tamasryszard.sztanka-toth@mdc-berlin.de, nikolaos.karaiskos@mdc-berlin.de

## Downstream stats

```{r parse_metadata, echo = F, val = F, eval = F}
metadata <- readStarLog(snakemake@input$star_log) %>%
    rbind(read_table2(snakemake@input$reads_type_out, col_names=c('observation', 'value'))) %>%
    # convert to millions
    mutate(value = round(value / 1e6, 2)) %>%
    {
        mutate(., input_reads = filter(., observation == 'input_reads')$value)
    } %>%
    mutate(label = ifelse(observation == 'input_reads', value, paste0(value, ' (', round(value/input_reads*100, 1), '%)'))) %>%
    dplyr::select(observation, label) %>%
    spread(observation, label) %>%
    dplyr::rename(as.utr = UTR,
                  intronic = INTRONIC,
                  intergenic = INTERGENIC,
                  as.cds = CODING) %>%
    # reorder columns
    dplyr::select(input_reads, uniq_mapped_reads, intergenic, intronic, as.cds, as.utr) %>%
    kable

metadata
```

## Saturation analysis

In order to know whether we would gain more from sequencing deeper, we downsampled the data (the final.bam file) to contain 10%, 20%... 90% reads, and then we created the DigitalExpression matrix (as in the normal dropseq pipeline).

This can give us insight, whether we have reached the saturation point (in terms of median umi per cell and median genes per cell) or whether we should sequence deeper.

Results of this are plotted below.

```{r read_summaries, echo = F}
downsampled_summaries <- snakemake@input[startsWith(names(snakemake@input), 'downsample')]

dge_data <- tibble(name = names(downsampled_summaries),
    filename = unlist(downsampled_summaries)) %>%
    separate(name, into = c(NA, 'run_mode', 'percentage'), sep='[.]') %>%
    mutate(content = map(filename, ~ read_csv(.) %>%
                         mutate(cell_bc = as.character(cell_bc)))) %>%
    unnest(content) %>%
    select(-filename) %>%
    rename(umis = total_counts,
        genes = n_genes_by_counts,
        reads = n_reads) %>%
    mutate(pcr = reads / umis)
```

```{r define_plotting_functions, echo =F }
cPalette = list('grey'= "#999999", 'light_orange'="#E69F00",
         'light_blue'="#56B4E9", 'green' = "#009E73",
         'yellow' = "#F0E442", 'blue'= "#0072B2", 'orange'="#D55E00",
         'pink'="#CC79A7")

clrs = c('umis'=cPalette$light_orange, 'pcr' = cPalette$light_blue, 'reads'=cPalette$green,
         'genes'=cPalette$pink, 'n_beads'= 'black')

median_clrs = clrs
names(median_clrs) = c(paste0('median_', names(clrs[-5])), names(clrs[5]))

plot_observations <- function(metric, run_mode_in, log_scale = T, ttl=''){
    pl <- dge_data %>%
        filter(run_mode == run_mode_in) %>%
        select(cell_bc, percentage, pcr, umis, reads, run_mode) %>%
        gather('obs', 'val', umis, pcr, reads) %>%
        filter(obs == metric) %>%
        mutate(percentage = as.integer(percentage)) %>%
        filter(percentage %in% c(20, 40, 60, 80, 100)) %>%
        group_by(percentage, obs, run_mode) %>%
        filter(between(val, quantile(val, 0.05, na.rm=T), quantile(val, 0.95, na.rm=T))) %>%
        mutate(percentage = factor(paste0(percentage, '%'), levels = c('20%', '40%', '60%', '80%', '100%'))) %>%
        ggplot(aes(val, fill = obs)) +
            geom_density() +
            facet_grid(percentage~run_mode, scales = 'free_y') +
            scale_fill_manual(values = clrs) +
            labs(x=ttl) +
            theme(legend.position='none', strip.background.x=element_blank(),
                  strip.text.x = element_blank(), strip.background.y = element_blank(),
                  strip.text.y = element_text(angle=0),
                  axis.ticks.y = element_blank(), axis.text.y = element_blank())
    if (log_scale){
        pl <- pl +
            scale_x_log10()+
            annotation_logticks(sides='b')
    }
    return(pl)
}

plot_deciled_data <- function(run_mode_in) {
    decile_dat <- dge_data %>%
        filter(run_mode == run_mode_in) %>%
        group_by(percentage) %>%
        mutate(cumsum_reads = cumsum(reads),
               decile_limit = sum(reads)/10,
               # put beads into deciles by number of reads
               decile = floor(cumsum_reads / decile_limit) + 1) %>%
        # get top 10 deciles, 11 is an artifact of rounding, last beads
        filter(decile < 11) %>%
        group_by(percentage, decile) %>%
        summarise(median_reads = median(reads),
                  median_genes = median(genes),
                  median_pcr = median(pcr),
                  median_umis = median(umis),
                  n_beads = n()) %>%
        gather('observation', 'value', median_reads:n_beads) %>%
        mutate(decile = factor(decile), 
               percentage = as.integer(percentage),
               observation = factor(observation, levels = c(
                          'median_reads', 
                          'median_umis', 
                          'median_genes',
                          'median_pcr', 'n_beads')))

    pl <- decile_dat %>%
        ggplot(aes(percentage, value, color= decile, fill = decile)) +
            geom_smooth(formula = y~log(x), size = 0.6, se=F) +
            geom_point(size=2,  color = 'black', pch=21) + 
            scale_x_continuous(breaks=seq(0, 100, 20)) +
            scale_fill_manual(values=parula(10)) + 
            scale_color_manual(values=parula(10)) + 
            facet_wrap(~observation, scales = 'free', ncol=2) +
            labs(y='', x='downsampling percentage') +
            guides(fill = guide_legend(nrow=3, byrow=T, override.aes = list(size=3))) +
            theme(legend.position = 'bottom', strip.background = element_blank())

    return(pl)
}

plot_data <- function(run_mode_in, obs_in, obs_name, umi_cutoff=c(1, 100, 200)){
    tibble(umi_cutoff = umi_cutoff,
                  dat = map(umi_cutoff, ~ filter(dge_data, umis > .))) %>%
        unnest(dat) %>%
        filter(run_mode == run_mode_in) %>%
        group_by(percentage, umi_cutoff) %>%
        summarise(median_reads = median(reads),
                  median_umis = median(umis),
                  median_genes = median(genes),
                  median_pcr = median(pcr),
                  n_beads = n()) %>%
        gather('observation', 'value', median_reads:n_beads) %>% 
        as_tibble() %>%
        mutate(percentage = as.integer(percentage),
               observation = factor(observation, levels = c(
                          'median_reads', 
                          'median_umis', 
                          'median_genes',
                          'median_pcr', 'n_beads'))) %>%
    mutate(umi_cutoff = factor(umi_cutoff)) %>%
    filter(observation == obs_in) %>%
    ggplot(aes(percentage, value, color = observation, fill=observation, linetype=umi_cutoff)) +
        scale_color_manual(values=median_clrs) +
        scale_fill_manual(values=median_clrs) +
        geom_smooth(formula = y~log(x), size = 0.6,se=F) +
        geom_point(size=2, color = 'black', pch=21) + 
        scale_x_continuous(breaks=seq(0, 100, 20), labels=paste0(seq(0, 100, 20), '%')) +
        labs(y=obs_name, x='downsampling percentage', color='', fill='', linetype='UMI cutoff') +
        theme(strip.background=element_blank(), legend.position='bottom',
              legend.key.width=unit(0.8, 'cm')) +
        guides(colour = 'none',fill='none',
               linetype = guide_legend(nrow=2, byrow=T, override.aes = list(size=1)))
}
```


## Histograms per run\_mode {.tabset}

```{r plot_histogram_of_observations, echo =F, fig.width=7, fig.height=4,results='asis'}
run_mode_names <- dge_data %$% 
    run_mode %>% unique()
    
for (run_mode_in in run_mode_names){
    umi_cutoff <- snakemake@params$run_modes[[run_mode_in]]$umi_cutoff
    cat(paste0('\n\n### ', run_mode_in, '\n\n'))
    print(plot_observations('umis', run_mode_in, ttl='# of UMIs per spatial unit'))
    print(plot_observations('reads', run_mode_in, ttl='# of reads per spatial unit'))
    print(plot_observations('pcr', run_mode_in, log_scale=F, ttl='reads / UMIs per spatial unit'))

}
```

## Median plots per run\_mode {.tabset}

```{r plot_median_values_of, echo =F, fig.width=7,fig.height=5,results='asis'}
run_mode_names <- dge_data %$% 
    run_mode %>% unique()
    
for (run_mode_in in run_mode_names){
    umi_cutoff <- snakemake@params$run_modes[[run_mode_in]]$umi_cutoff
    cat(paste0('\n\n### ', run_mode_in, '\n\n'))

    print(plot_data(run_mode_in, obs_in = 'median_reads',
                    obs_name = 'median reads\nper spatial unit', umi_cutoff = c(1, umi_cutoff)))
    print(plot_data(run_mode_in, obs_in = 'median_umis',
                    obs_name = 'median UMIs\nper spatial unit', umi_cutoff = c(1, umi_cutoff)))
    print(plot_data(run_mode_in, obs_in = 'median_pcr',
                    obs_name = 'median reads/UMIs\nper spatial unit', umi_cutoff = c(1, umi_cutoff)))
}
```

## Deciled median plots per run\_mode {.tabset}

```{r plot_deciled_median_values_of, echo =F, fig.width=7,fig.height=7,results='asis'}
run_mode_names <- dge_data %$% 
    run_mode %>% unique()
    
for (run_mode_in in run_mode_names){
    umi_cutoff <- snakemake@params$run_modes[[run_mode_in]]$u5i_cutoff
    cat(paste0('\n\n### ', run_mode_in, '\n\n'))

    print(plot_deciled_data(run_mode_in))
}
```
---
header-includes:
    - \usepackage{float}
    - \usepackage[table]{xcolor}
output:
    html_document:
        toc: true
        toc_depth: 6
classoption: landscape
geometry: margin=0.5cm
version: 0.1.1
author: Tamas Ryszard Sztanka-Toth, Nikolaos Karaiskos
email: tamasryszard.sztanka-toth@mdc-berlin.de, nikolaos.karaiskos@mdc.berlin.de
license: GPL
title: Sample overview
pagetitle: Sample overview
date: "`r format(Sys.time(),'%d/%m/%y')`"
---

```{r knitr_options, include=FALSE, cache=FALSE}
knitr::opts_chunk$set(
  cache = F,
  autodep = TRUE,
  message = FALSE,
  warning = FALSE,
  comment = NA
)

options(knitr.table.format ='markdown')
```

```{r functions, echo = F}
readStarLog <- function(log_file){

		out = list()
		lines = readLines(log_file)
	
		out$input_reads = (lines[6] %>% strsplit('\t') %>% unlist)[2] %>% as.integer

		out$uniq_mapped_reads = (lines[9] %>% strsplit('\t') %>% unlist)[2] %>% as.integer

		#out$avg_length = (lines[11] %>% strsplit('\t') %>% unlist)[2] %>% as.numeric
		
        tibble(observation=names(out), value=unlist(unname(out)))
	}
```

```{r load_projects_puck_info, echo=F}
library(tidyverse)
library(magrittr)
metadata <- read_csv(snakemake@input[[1]])
```

```{r collect_data, echo = F}
root_dir <- snakemake@config$root_dir
read_metrics <- metadata %>%
    select(project_id, sample_id, puck_id, species, sequencing_date) %>%
    mutate(star_log = paste0(root_dir, '/projects/', project_id, '/processed_data/', sample_id, '/illumina/complete_data/star_Log.final.out'),
           read_types =paste0(root_dir,'/projects/', project_id, '/processed_data/', sample_id, '/illumina/complete_data/split_reads/read_type_num.txt')) %>%
    mutate(star_log = map(star_log,
                          ~ readStarLog(.))) %>%
    unnest(star_log) %>%
    mutate(read_types = map(read_types,
                            ~ read_table2(., col_names=c('rt_obs', 'rt_value')))) %>%
    unnest(read_types) %>%
    mutate(rt_obs = tolower(rt_obs)) %>%
    spread(rt_obs, rt_value) %>%
    spread(observation, value)
```

```{r show_sample_table, echo = F}
library(kableExtra)
to_table <- read_metrics %>%
    mutate(um_r = uniq_mapped_reads) %>%
    gather('obs', 'val', intergenic, amb, coding, intronic, utr) %>%
    mutate(val_p = round(val / um_r, 2),
           val = round(val / 1e6, 2),
           # add ratio in paranthesis if obs is not cds
           val = paste0(val, ' (', val_p, ')'),
           uniq_mapped_reads = round(uniq_mapped_reads / 1e6, 2),
           input_reads = round(input_reads / 1e6, 2),
           uniq_mapped_reads = paste0(uniq_mapped_reads, ' (', round(uniq_mapped_reads / input_reads, 2), ')')) %>%
    select(-um_r, -val_p) %>%
    spread(obs, val) %>%
    arrange(species) %>%
    select(sample_id, puck_id, species, sequencing_date, input_reads, uniq_mapped_reads, coding, utr, intergenic, intronic, amb) %>%
    rename(uniq_m = uniq_mapped_reads,
           input_r = input_reads,
           cds = coding)
```

```{r load_strand_info, echo = F}
strand_info <- metadata %>%
    select(project_id, sample_id, puck_id, species, sequencing_date) %>%
    mutate(filename = paste0(root_dir, '/projects/', project_id, '/processed_data/', sample_id, '/illumina/complete_data/split_reads/strand_type_num.txt'),
           content = map(filename, ~read_table2(., col_names = c('obs', 'num')))) %>%
    unnest(content) %>%
    select(-filename, project_id) %>%
    group_by(sample_id) %>%
    mutate(num_sum = sum(num),
           num_ratio = round(num / num_sum, 2),
           num = round(num / 1e6, 2),
           num = paste0(num, ' (', num_ratio, ')')) %>%
    select(-num_ratio, -num_sum) %>%
    spread(obs, num)
```

```{r load_barcode_metadata, echo = F}
umi_cutoffs <- c(1, 10, 50, 100)

load_filter_dge <- function(x, y){
    read_table2(x, skip=6) %>%
        filter(NUM_TRANSCRIPTS > y)
}

read_dge_summary <- function(filename){
    tibble(umi_cutoff = umi_cutoffs, filename=filename) %>%
        mutate(dat = map2(filename, umi_cutoff, load_filter_dge)) %>%
        select(-filename) %>%
        unnest(dat) %>%
        group_by(umi_cutoff) %>%
        summarise(
                median_umi = median(NUM_TRANSCRIPTS),
                median_reads = median(NUM_GENIC_READS),
                median_genes = median(NUM_GENES),
                median_pcr = median(round(NUM_GENIC_READS / NUM_TRANSCRIPTS, 1)),
                mean_umi = as.integer(mean(NUM_TRANSCRIPTS)),
                mean_reads = as.integer(mean(NUM_GENIC_READS)),
                mean_genes = as.integer(mean(NUM_GENES)),
                num_beads = n())
        
}

barcode_metadata <- metadata %>%
    select(project_id, sample_id, puck_id, species, sequencing_date) %>%
    mutate(filename = paste0(root_dir, '/projects/', project_id, '/processed_data/', sample_id, '/illumina/complete_data/dge/')) %>%
    mutate(filename = ifelse(file.exists(paste0(filename, 'dge_all_summary.txt')),
                             paste0(filename, 'dge_all_summary.txt'),
                             paste0(filename, 'dge_all_cleaned_summary.txt')),
           content = map(filename, ~read_dge_summary(.))) %>%
    select(-filename, -project_id) %>%
    unnest(content)
```

## Overview

We show here downstream metadata for each experiment performed in the sts project. There are three types of tables:

* Read information table: containing the parsed output of mapping, such as input read number, uniquely mapped read number etc.
* Expression summary table: containing median number of umis, genes, reads (and mean) per bead for each sample. This is done after applying a UMI filter of 1, 10, 50, 100.
* Strand information table: containing the numbers for reads mapping to the correct strand

Each table has the following 4 columns: sample\_id, puck\_id, species, sequencing\_date

### Table column description

__Read information table__

* input\_r: number of input reads (millions) from the flowcell
* uniq\_m: number of uniquely mapped reads (millions). In parantheses ratio to input\_r
* cds, utr, intergenic, intronic, amb: coding, utr, intergenic, intronic and ambient (overlapping genes on both strands, or cannot be assigned to a single gene part). In millions, in parantheses ratio to uniq\_m.

__Expression summary tables__

All columns here are in raw counts. We have mean and median for UMIs, genes, reads (all per bead). Median pcr is the median of reads/umi (per bead).

__Strand information table__

Here there are 6 columns: minus\_AMB,  minus\_minus, minus\_plus, plus\_AMB, plus\_minus, plus\_plus. The first part is the position of the read (plus or minus strand) the second is the position of the mapped gene. AMB means that the mapped gene is ambient (overlapping genes on different strand) or that the read is intergenic.

## Tables by species containing sequencing metadata


```{r print_by_species, echo = F, results = 'asis'}
for(s in unique(to_table$species)){
    cat(paste0('### ', s, ' samples'))
    cat('\n')

    cat('#### Read information table\n') 
    to_table %>%
        filter(species == s) %>%
        kable("html") %>%
        kable_styling('striped', font_size=12)  %>%
        row_spec(row=0, bold=T) %>%
        print

    cat('\n')
    cat('[Back to top](#)\n\n')

    cat('#### Expression summary tables\n')

    for(cutoff in umi_cutoffs){
        cat(paste0('##### UMI cutoff: ', cutoff))
        cat('\n')

        barcode_metadata %>%
            filter(species == s, umi_cutoff == cutoff) %>%
            kable("html") %>%
            kable_styling('striped', font_size=12)  %>%
            row_spec(row=0, bold=T) %>%
            print

        cat('\n')
        cat('[Back to top](#)\n\n')

    }

    cat('#### Strand information table\n')

    strand_info %>%
        filter(species == s) %>%
        kable("html") %>%
        kable_styling('striped', font_size=12)  %>%
        row_spec(row=0, bold=T) %>%
        print

    cat('\n')
    cat('[Back to top](#)\n\n')
}
```


---
output:
  html_document:
    toc: true
    self_contained: yes
    toc_float: true
    toc_depth: 4
    theme: flatly
    highlight: tango
version: 0.3.0
author: Tamas Ryszard Sztanka-Toth, Nikolaos Karaiskos
email: tamasryszard.sztanka-toth@mdc-berlin.de, nikolaos.karaiskos@mdc.berlin.de
license: GPL
---

```{r knitr_options, include=FALSE, cache=FALSE}
knitr::opts_chunk$set(
  autodep = TRUE,
  message = FALSE,
  warning = FALSE,
  comment = NA
)
```

```{r libraries, include = F, cache=F}
library(tidyverse)
library(yaml)
library(knitr)
library(magrittr)
library(kableExtra)
library(cowplot)

theme_set(theme_cowplot(18))

cpalette <- list('orange' = '#D55E00', 'blue' = '#0072B2', 'green' = '#009E73', 'black' = '#000000', 'yellow' = '#F0E442', 
				 'grey' = '#999999', 'light_orange' = "#E69F00", 'light_blue' = "#56B4E9")

readStarLog <- function(log_file){

		out = list()
		lines = readLines(log_file)
	
		out$input_reads = (lines[6] %>% strsplit('\t') %>% unlist)[2] %>% as.integer

		out$uniq_mapped_reads = (lines[9] %>% strsplit('\t') %>% unlist)[2] %>% as.integer

        out$multi_mapped_reads = (lines[24] %>% strsplit('\t') %>% unlist)[2] %>% as.integer

		out$avg_mapped_length = (lines[11] %>% strsplit('\t') %>% unlist)[2] %>% as.numeric

        out$unmapped_too_short = (lines[31] %>% strsplit('\t') %>% unlist)[2] %>% as.integer
		
        tibble(observation=names(out), value=unlist(unname(out)))
	}
```

```{r read_sample_information_tables, echo = F}
run_mode_settings <- bind_rows(snakemake@params$run_mode) %>%
    nest(umi_cutoff = c(umi_cutoff)) %>%
    bind_cols(tibble(run_mode = names(snakemake@params$run_mode))) %>%
    gather('variable_name', 'value', -run_mode) %>%
    spread(run_mode, value)

run_mode_names <- run_mode_settings[-1] %>% names
run_mode_n <- length(run_mode_names)
normal_run_mode_names <- run_mode_settings %>%
    gather('run_mode','val', -variable_name) %>%
    filter(variable_name == 'mesh_data') %>%
    unnest(val) %>%
    filter(val == 0) %$%
    run_mode

meshed_run_mode_names <- run_mode_settings %>%
    gather('run_mode','val', -variable_name) %>%
    filter(variable_name == 'mesh_data') %>%
    unnest(val) %>%
    filter(val == 1) %$%
    run_mode

sample_info <- bind_rows(snakemake@params$sample_info) %>%
    select(species, sequencing_date, investigator, experiment) %>%
    unique %>%
    gather('info_name', 'value') %>%
    add_row(info_name = 'project_id',
            value = snakemake@wildcards$project, .before=1)%>%
    add_row(info_name = 'sample_id',
            value = snakemake@wildcards$sample, .after=1)%>%
    add_row(info_name = 'puck_id',
            value = snakemake@wildcards$puck, .after=2)

cPalette = list('grey'= "#999999", 'light_orange'="#E69F00",
         'light_blue'="#56B4E9", 'green' = "#009E73",
         'yellow' = "#F0E442", 'blue'= "#0072B2", 'orange'="#D55E00",
         'pink'="#CC79A7")

clrs = c('umis'=cPalette$light_orange, 'genes' = cPalette$light_blue, 'reads'=cPalette$green,
         'pcr'=cPalette$pink)

summarised_clrs = c(clrs, 'black')
names(summarised_clrs) = c(paste0('median_', names(clrs)), 'n_beads')

nucl_clrs = c('A'='#F5C900',
              'C'='#F55D59',
              'T'='#3AA861',
              'G'='#7772F5',
              'N'='#999999')

# loading the dge summaries
dge_summary <- tibble(run_mode = run_mode_names,
                      obs_df = snakemake@input[paste0(run_mode, '.dge_summary')]) %>%
    unnest(obs_df) %>%
    mutate(dat = map(obs_df, ~ read_csv(.) %>%
        mutate(cell_bc = as.character(cell_bc)))) %>%
    unnest(dat) %>%
    select(-obs_df) %>%
    rename(pcr = reads_per_counts,
           umis = total_counts,
           genes = n_genes_by_counts,
           reads = n_reads) %>%
    group_by(run_mode) %>%
    filter(grepl('^[0-9]+|^[ACTGAN]+$', cell_bc))
```

## Overview {.tabset}

### Sample info

```{r show_sample_info, echo = F}
sample_info %>%
    kbl(col.names=NULL) %>%
    kable_classic_2(full_width=F, position='left')
```

### Run modes

This sample was processed using the following run modes, and run mode variables

```{r plot_run_modes, echo = F}
col_names_to_show <- colnames(run_mode_settings)[-1]

run_mode_settings %>%
    kbl(col.names = c('', col_names_to_show)) %>%
    kable_classic_2(full_width=F, position='left') %>%
    add_header_above(c(" ", "Run modes" = length(col_names_to_show))) %>%
    column_spec(1, border_right=T) %>%
    row_spec(0, bold=T)
```

```{r load_read_statistics, echo =F}
rRNA_stats <- read_table2(snakemake@input$ribo_log, col_names=c('observation', 'value')) %>%
    spread(observation, value) %>%
    mutate(mapped_to_rRNA = aligned_reads) %>%
    gather('observation', 'value') %>%
    filter(observation == 'mapped_to_rRNA') %>%
    mutate(value = ifelse(value == 'None', 0, value))
 
read_stats <- tibble(star_logs = snakemake@input$star_log,
                     reads_type_out = snakemake@input$reads_type_out) %>%
    mutate(trimmed = ifelse(grepl('polyA_adapter_trimmed', reads_type_out, fixed=TRUE),
                            'polyA_adapter_trimmed', 'untrimmed')) %>%
    mutate(star_log_dat = map(star_logs, ~ readStarLog(.)),
           reads_type_out_dat = map(reads_type_out, ~ read_table2(.,
                                                                  col_names = c('observation',
                                                                                'value')))) %>%
    mutate(rRNA_stats_dat = map(star_logs, ~ rbind(rRNA_stats))) %>%
    select(-star_logs, -reads_type_out) %>%
    gather('dat_type', 'dat', -trimmed) %>%
    select(-dat_type) %>%
    unnest(dat) %>%
    group_by(trimmed) %>%
    #rbind(rRNA_stats) %>%
    # convert to millions
    mutate(value = ifelse(observation != 'avg_mapped_length',
                          round(value / 1e6, 2), value)) %>%
    # add input_reads as a column per group
    spread(observation, value) %>%
    mutate(inp_reads = input_reads) %>%
    gather('observation', 'value', -trimmed, -inp_reads) %>%
    mutate(label = ifelse(observation %in% c('input_reads', 'avg_mapped_length'),
                          value, paste0(value, ' (', round(value/inp_reads*100, 1), '%)'))) %>%
    dplyr::select(observation, label) %>%
    spread(observation, label) %>%
    dplyr::rename(as.utr = UTR,
                  intronic = INTRONIC,
                  intergenic = INTERGENIC,
                  ambiguous = AMB,
                  as.cds = CODING) %>%
    # reorder columns
    dplyr::select(input_reads, uniq_mapped_reads, avg_mapped_length, multi_mapped_reads, 
                  unmapped_too_short, intergenic, intronic, as.cds, ambiguous, as.utr,
                  mapped_to_rRNA) %>%
    gather('metric', 'value', -trimmed) %>%
    spread(trimmed, value) %>%
    # change order
    slice(5, 10, 2,3,7,6,1,4,9,11,8)

map_types = colnames(select(read_stats, -metric))
```

### Mapping statistics

The sample was mapped using `r map_types` reads. The mapping statistics are shown here for each method:

```{r show_read_stats, echo = F}
col_names_to_show <- colnames(read_stats)[-1]
read_stats %>%
    kbl(col.names = c('', col_names_to_show)) %>%
    kable_classic_2(full_width=F, position='left') %>%
    add_header_above(c(" ", "Mapping mode" = length(col_names_to_show))) %>%
    column_spec(1, border_right=T) %>%
    row_spec(0, bold=T) %>%
    add_indent(c(3,4,5,6,7)) %>%
    footnote(general ='All values, except avg_mapped_length, shown in millions')
```

### Summarised metrics over beads

```{r calculate_summarised_metrics_over_beads, echo = F}
col_names_to_show <- colnames(run_mode_settings)[-1]

dge_summary %>%
    group_by(run_mode) %>%
    summarise(sum_reads = sum(reads, na.rm=T),
              median_reads = median(reads, na.rm=T),
              median_umis = median(umis, na.rm=T),
              median_pcr = median(pcr, na.rm=T),
              median_genes = median(genes, na.rm=T),
              n_beads = n()) %>%
    mutate(sum_reads = paste0(round(sum_reads / 1e6, 2), ' (1e6)')) %>%
    gather('obs', 'value', -run_mode) %>%
    spread(run_mode, value) %>%
    kbl(col.names = c('', col_names_to_show)) %>%
    kable_classic_2(full_width=F, position='left') %>%
    add_header_above(c(" ", "Run modes" = length(col_names_to_show))) %>%
    column_spec(1, border_right=T) %>%
    row_spec(0, bold=T)

```

## QC plots

Each of the QC plots we show on a per run mode basis, to see if there are any downstream differences based on the run mode variable settings.

### 'Knee'-plot {.tabset}

Below we plot a so called 'Knee-plot': on the y-axis is the Cummulative sum of reads, on the x-axis are the bead barcodes sorted by number of reads. For single-cell samples, this plot tells you roughly how many beads are in the sample.

```{r knee_plot, echo = F, fig.height=4, fig.width=7, results='asis'}
read_counts <- dge_summary %>%
    select(run_mode, cell_bc, reads) %>%
    mutate(reads_cumsum = cumsum(reads),
           ix = 1:n())

for (run_mode_in in unique(dge_summary$run_mode)){
    cat (paste0('\n\n #### ', run_mode_in, '\n\n'))

    pl_knee <- read_counts %>%
        filter(run_mode == run_mode_in) %>%
        ggplot(aes(ix, reads_cumsum)) +
            geom_line() +
            labs(x='Beads sorted by number of reads', y='Cummulative\nsum of reads')

    print(pl_knee)
    cat('\n\n')
}
```

### Umi-cutoff plots

```{r create_umi_cutoffs, echo = F, fig.width = 10, fig.height=5}
umi_cutoffs <- seq(10, 20000, 10)

dge_summary <- dge_summary %>%
    group_by(run_mode) %>%
    mutate(reads_cumsum = cumsum(reads)) %>%
    mutate(quartile = cut(reads_cumsum,
                          breaks = 4,
                          include.lowest= T,
                          labels = c('Q1', 'Q2', 'Q3', 'Q4'))) %>%
    select(-reads_cumsum)

summarise_dge_summary <- function(umi_cutoff){
    dge_summary %>%
        filter(umis > umi_cutoff) %>%
        summarise(median_reads = median(reads),
                  median_umis = median(umis),
                  median_genes = median(genes),
                  median_pcr = median(pcr),
                  n_beads = n())
}

umi_cutoff_data <- tibble(umi_cutoffs = umi_cutoffs) %>%
    mutate(dat = map(umi_cutoffs, ~ summarise_dge_summary(.))) %>%
    unnest(dat)

```


```{r plot_umi_cutoff_plot, echo = F, fig.height=7.5, fig.width=12}
metric_names <- list('umis' = 'UMIs', 'genes' = 'genes', 'pcr' = 'reads / UMIs',
                     'reads'='reads')
umi_cutoff_plot <- function(metric, y_log_transform = F, legend_pos='none'){
    y_label <- strsplit(metric, '_')
    y_label <- paste0(y_label[[1]][1], ' ', metric_names[[y_label[[1]][2]]])

    pl <- umi_cutoff_data %>%
        gather('obs', 'value', -umi_cutoffs, -run_mode) %>%
        filter(obs == metric) %>%
        ggplot(aes(umi_cutoffs, value, color=obs, linetype = run_mode)) +
            geom_line() +
            scale_color_manual(values=summarised_clrs) +
            scale_x_log10(breaks = c(1e1, 1e2, 1e3, 1e4)) +
            theme(strip.background=element_blank(), strip.text.x=element_blank(),
                  text = element_text(size=18, face = 'plain'), legend.position=legend_pos,
                  legend.title = element_blank(), legend.margin = margin(t=0.5, b=0.5,l=0.5,
                                                                         r=0.5, unit = 'cm'),
                  legend.spacing = unit(0, 'cm'), plot.margin = unit(c(1,1,1,1), "lines")) +
            labs(color='', linetype='', y=paste0(y_label, '\nper spatial unit'), x='minimum UMI') +
            guides(linetype = guide_legend(override.aes = list(size = 1), ncol=2),
                   color = guide_legend(override.aes = list(size = 2)))

    if (y_log_transform){
        pl <- pl + scale_y_log10() +
            annotation_logticks(sides='bl')
    } else {
        pl <- pl +
            annotation_logticks(sides='b')
    }

    return(pl)
}

pl1 <- umi_cutoff_plot('n_beads', y_log_transform = T)
pl2 <- umi_cutoff_plot('median_reads')
pl3 <- umi_cutoff_plot('median_genes')
pl4 <- umi_cutoff_plot('median_umis')
pl5 <- umi_cutoff_plot('median_pcr', legend_pos = 'right')

plt_legend <- get_legend(pl5)

plot_grid(pl1, pl2, pl3,
        pl4, pl5 + theme(legend.position='none'), plt_legend,
        align='vh', hjust=-1, labels="", nrow=3)
```

### Histogram of metrics over beads {.tabset}

Next we show mertics such as number of UMIs, genes, reads and pcr per physical spot. We further distinguish between each run mode, showing one histogram for each. 

```{r plot_n_reads_bead_hist, echo = F, fig.width=10, fig.height=4.5, dpi=300, results='asis'}
library(scales)
to_plot <- dge_summary %>%
    gather('obs', 'val', -cell_bc, -quartile, -run_mode)

min_difference <- to_plot %>% group_by(obs) %>%
    filter(!is.na(val)) %>%
    summarise(min_val = min(val,na.rm=T ),
              max_val = max(val,na.rm=T),
              difference = max_val - min_val) %>%
    summarise(min_difference = min(difference)) %$%
    min_difference

metric_plot <- function(run_mode_in, metric, legend_pos='none'){
    y_label <- metric_names[[metric]]
    to_plot <- to_plot %>%
        filter(obs == metric)

    pl <-  to_plot %>%
       filter(run_mode == run_mode_in) %>% 
       ggplot(aes(x = val, fill=obs)) +
            geom_histogram(bins=100) +
            scale_x_log10(breaks = c(1, 1e1, 1e2, 1e3, 1e4, 1e5, 1e6),
                          labels=trans_format('log10',math_format(10^.x))) +
            scale_fill_manual(values=clrs) +
            labs(fill='', x=paste0('# of ', y_label),
                 y=paste0('# of\n',
                          'spatial units')) +
            theme(legend.position = legend_pos, strip.background=element_rect(fill='white'),
                  text = element_text(size =18 ))

    if(min_difference > 1){
        pl <- pl + 
            annotation_logticks(sides='b')
    }

    return(pl)
}
for (run_mode in unique(dge_summary$run_mode)){
    cat (paste0('\n\n #### ', run_mode, '\n\n'))

    pl1 = metric_plot(run_mode, 'reads')
    pl2 = metric_plot(run_mode, 'pcr')
    pl3 = metric_plot(run_mode, 'genes')
    pl4 = metric_plot(run_mode, 'umis')

    print(plot_grid(pl1, pl2, pl3,
            pl4, align='vh', hjust=-1, labels="", nrow=2))

    cat('\n\n')
}
```

### Nucleotide distribution per beads {.tabset}

Next we bin the data based on reads into quartile. For each run\_mode the data is divided into 4 beads, by reads. This means, that the first bin will contain beads which account 25% of the reads, the second will contain beads which account for the second 25% of reads and so on. 

For each run mode we plot the nucleotide distribution per quartile.

**You can pick the run\_mode to be shown from the dropdown list**

**Only not-meshed run\_mode(s) are shown**

```{r plot_nucl_freq, results='asis', fig.width=10, fig.height = 4.5, echo=F, dpi=300}
plot_nucl_freq <- function(dge_s){
    cell_bc_len <- nchar((dge_s %$% cell_bc)[1])

    nucls <- dge_s %$%
        cell_bc %>% strsplit("")

    nucls <- dge_s %>%
        select(cell_bc, quartile, run_mode) %>%
        bind_cols(tibble(nucl=nucls)) %>%
        unnest(nucl) %>%
        # group by run_mode and cell_bc, so each nt is on the correct pos
        group_by(run_mode, cell_bc) %>%
        mutate(pos = paste0('pos_', 100+1:n())) %>%
        group_by(run_mode, pos, nucl, quartile) %>%
        summarise(nucl_count = n()) %>%
        ungroup() %>%
        tidyr::complete(run_mode, pos, nucl, quartile, fill=list(nucl_count=0))

    lbl_df <- dge_s %>%
        ungroup() %>%
        group_by(quartile) %>%
        summarise(lbl = n()) %>%
        mutate(lbl = paste0(quartile, ' (n=', lbl, ')'))

    lbls <- lbl_df$lbl
    names(lbls) <- lbl_df$quartile

    nucls %>%
        ggplot(aes(pos, nucl_count, fill = nucl)) +
            geom_bar(stat='identity', position='dodge') +
            scale_fill_manual(values=nucl_clrs) +
            scale_x_discrete(labels=seq(1, cell_bc_len, 1)) +
            facet_wrap(~quartile, ncol=2, scales='free_y',
                       labeller=as_labeller(lbls)) +
            labs(x='nucleotide position in the barcode',
                 y='nucleotide count')+
            theme(text = element_text(size =18 ))
}

for (run_mode_in in unique(normal_run_mode_names)){
    cat (paste0('\n\n #### ', run_mode_in, '\n\n'))
    dge_summary %>%
        filter(run_mode == run_mode_in) %>%
        plot_nucl_freq() %>%
        print()
}
```

### Shannon entropy and string compression {.tabset}

```{r plot_shannon, echo = F, fig.height=5, fig.width=10, results='asis', dpi=300}
plot_shannon_scompression <- function(run_mode_in){
    dat <- dge_summary %>%
        filter(run_mode == run_mode_in) %>%
        select(cell_bc, run_mode, ends_with('compression') | ends_with('entropy')) %>%
        gather('observation', 'value', -cell_bc, -run_mode) %>%
        # replace first _ with |, so we can separate later
        mutate(observation = str_replace(observation, '_', ' '),
               observation = str_replace(observation, '\\s', '|'),
               observation = str_replace(observation, ' ', '_')) %>%
        separate(observation, into = c('type', 'observation'), sep = '\\|')

    pl1 <- dat %>%
        filter(observation == 'entropy') %>%
        ggplot(aes(value, fill = type)) +
            geom_histogram(bins=30, color='black', position='dodge') +
            scale_fill_manual(values=c(cpalette$grey, cpalette$orange),labels=c('theoretical', 'observed')) +
            labs(fill='', y='# of barcodes', x='Shannon entropy of barcodes') +
            theme(legend.position = c(0.1, 0.8), text = element_text(size=18, face='plain'))

    pl2 <- dat %>%
        filter(observation == 'compression') %>%
        ggplot(aes(value, fill = type)) +
            geom_histogram(bins=30, color='black', position='dodge') +
            scale_fill_manual(values=c(cpalette$grey, cpalette$orange)) +
            scale_x_continuous(limits=c(0, NA)) +
            labs(fill='', y='# of barcodes', x='Length of barcodes after compression') +
            theme(text = element_text(size=18, face='plain'), legend.position='none')

    plt_legend <- get_legend(pl2)

    plot_grid(pl1, pl2,
              align='vh', hjust=-1, labels="", nrow=2)

}

for (run_mode in unique(dge_summary$run_mode)){
    cat (paste0('\n\n #### ', run_mode, '\n\n'))

    print(plot_shannon_scompression(run_mode))
}
```

```{r prepare_spatial, echo =F}
# barcode file attached at the python level
is_spatial <- snakemake@params$is_spatial
```

```{r spatial_header, eval=is_spatial, echo=F, results='asis'}
cat('## Spatial QC {.tabset}\n\n')
```

```{r prepare_spatial_data, eval=is_spatial, echo =F}
run_mode_info <- dge_summary %>%
    select(run_mode) %>%
    summarise(n_cells=n()) %>%
    mutate(def_plot_bead_size = ifelse(n_cells > 5000, 0.4, 0.75)) %>%
    mutate(def_plot_bead_size = ifelse(n_cells > 10000, 0.1, def_plot_bead_size)) %>%
    mutate(def_plot_bead_size = ifelse(n_cells > 25000, 0.05, def_plot_bead_size))

plot_bead_size <- run_mode_settings %>%
    filter(variable_name %in% c('mesh_data', 'mesh_spot_diameter_um')) %>%
    gather('run_mode', 'variable_value', -variable_name) %>%
    unnest(variable_value) %>%
    add_row(variable_name = 'spot_diameter_um', run_mode = run_mode_names,
            variable_value = snakemake@params$puck_variables$spot_diameter_um) %>%
    spread(variable_name, variable_value) %>%
    inner_join(run_mode_info, by='run_mode') %>%
    mutate(mesh_data = as.logical(mesh_data)) %>%
    group_by(run_mode) %>%
    mutate(plot_bead_size = min(ifelse(mesh_data, mesh_spot_diameter_um / 40,
                              spot_diameter_um / 40), def_plot_bead_size)) %>%
    select(run_mode, plot_bead_size)

x_limits <- dge_summary$x_pos %>% {c(min(.), max(.))}
y_limits <- dge_summary$y_pos %>% {c(min(.), max(.))}


ratio <- (x_limits[2] - x_limits[1] ) / (y_limits[2] - y_limits[1])

puck_width_um <- snakemake@params$puck_variables$width_um
scale_factor <- ifelse(puck_width_um < 3000, 2, 3)
mm_dist <- max(10^scale_factor, round(puck_width_um/3, digits =-scale_factor))
mm_diff <- mm_dist / 1000

x_mm_breaks <- seq(0, puck_width_um, mm_dist)
x_mm_breaks <- paste0(x_mm_breaks * mm_diff / mm_dist, 'mm')
y_mm_breaks <- seq(0, puck_width_um / ratio, mm_dist)
y_mm_breaks <- paste0(y_mm_breaks * mm_diff / mm_dist, 'mm')

px_to_um <- (x_limits[2] - x_limits[1]) / snakemake@params$puck_variables$width_um

x_breaks <- seq(x_limits[1], x_limits[2], px_to_um * mm_dist)
y_breaks <- seq(y_limits[1], y_limits[2], px_to_um * mm_dist)

plot_spatial_qc <- function(variable, dat_in, plot_bead_size=0.1, ttl='', cut_dist=T, min_value=0, is_hexagonal=F){
    spatial_qc_dat <- dat_in %>%
        mutate(qnt = quantile(val, 0.9, na.rm=T)) %>%
        arrange(val)

    if(cut_dist){
    spatial_qc_dat <- spatial_qc_dat %>%
        mutate(val = ifelse(val > qnt, qnt, val))
    }

    limits <- c(min_value, quantile(spatial_qc_dat$val, 0.91, na.rm=T))

    if(is_hexagonal){
        pl <- spatial_qc_dat %>%
            ggplot(aes(x_pos, y_pos, fill = val, group=1)) +
                geom_hex(stat='identity', color='gray95') +
                guides(fill = guide_colorbar(barheight = 10)) + 
                labs(fill=ttl, x='', y='')
        if(cut_dist){
            pl <- pl +
                scale_fill_viridis_c(option =  "magma", limits = limits)
        } else {
            pl <- pl + 
                scale_fill_viridis_c(option =  "magma") 
        }
    } else{
        pl <- spatial_qc_dat %>%
            ggplot(aes(x_pos, y_pos, color = val)) +
                geom_point(size=plot_bead_size) + 
                guides(color = guide_colorbar(barheight = 10)) + 
                labs(color=ttl, x='', y='')
        if(cut_dist){
            pl <- pl +
                scale_color_viridis_c(option =  "magma", limits = limits)
        } else {
            pl <- pl + 
                scale_color_viridis_c(option =  "magma") 
        }
    }

    pl <- pl + 
        scale_x_continuous(labels = x_mm_breaks, breaks = x_breaks, limits=x_limits) +
        scale_y_continuous(labels = y_mm_breaks, breaks = y_breaks, limits=y_limits) +
        coord_fixed()+
        theme(panel.background = element_rect(fill = 'gray95'), 
             legend.spacing = unit(0.1, 'cm'),
            axis.line = element_line(color = 'black'),
            text = element_text(color='black', size=18))

    return(pl)
}

qc_vars <- c('n_joined', 'reads', 'genes', 'umis', 'pcr',
               'pct_counts_mt', 'exact_entropy',
               'exact_compression')

qc_var_info <- dge_summary %>%
    ungroup() %>%
    select(-cell_bc, -x_pos, -y_pos, -total_counts_mt, -run_mode,
           -theoretical_entropy, -theoretical_compression, -quartile) %>%
    unique() %>%
    colnames() %>%
    tibble(qc_vars = .) %>%
    inner_join(tibble(qc_vars = qc_vars,
               ttl = c('# beads joined\nper spatial unit',
                       '# of reads\nper spatial unit',
                       '# of genes\nper spatial unit',
                       '# of UMIs\nper spatial unit',
                       'reads / UMIs\nper spatial unit',
                       '% mt counts\nper spatial unit',
                       'Shannon\nentropy\nper spatial unit',
                       'barcode length\nafter\ncompression\nper spatial unit'),
               cut_dist = c(F, T, T, T, T, T, F, F),
               min_value = c(0, 0, 0, 0, 1, 0, 0, 0)),
               by='qc_vars') %>%
    mutate(qc_vars = factor(qc_vars, levels= qc_vars))
```


```{r plot_per_run_mode, echo =F, results='asis', fig.width=7, fig.height=6, eval=is_spatial, dpi=300}
for (run_mode_name in run_mode_names){
    cat(paste0('\n\n### ', run_mode_name, '\n\n'))

    for (qc_var in qc_var_info$qc_vars){
        pbs = filter(plot_bead_size, run_mode == run_mode_name)$plot_bead_size
        dat_in = dge_summary %>%
            filter(run_mode == run_mode_name) %>%
            ungroup() %>%
            select(x_pos, y_pos, matches(qc_var)) %>%
            rename(val = qc_var) %>%
            filter(!is.na(val))

        if (qc_var %in% colnames(dge_summary) & nrow(dat_in) > 0){
            qc_var_info %>%
                filter(qc_vars == qc_var) %>%
                {
                    plot_spatial_qc(
                       qc_var,
                       dat_in,
                       plot_bead_size = pbs,
                       ttl = .$ttl,
                       cut_dist = .$cut_dist,
                       min_value = .$min_value,
                       is_hexagonal=snakemake@params$run_modes[[run_mode_name]]$mesh_type == 'hexagon')
                } %>%
                print()
        }

    }
}
```

Initialization
==============

Initializing using required arguments
-------------------------------------

.. include:: shared/spacemake_init.rst

Optional arguments
------------------

The `spacemake init` command takes the following optional arguments:

``root_dir``
    The ``root_dir`` for the spacemake instance. Defaults to ``.``, the directory in which `spacemake init` is ran.

``temp_dir``
    Path to the temporary directory, defaults to ``/tmp``.

``download_species``
    If set, spacemake will download the genome (.fa) and annotation (.gtf) files for mouse and human (from gencode, as specified `here <https://github.com/rajewsky-lab/spacemake/blob/master/spacemake/data/config/species_data_url.yaml>`_.

Hence, the complete `spacemake init` command looks like this::
    
    spacemake init \
      --root_dir ROOT_DIR \             # optional
      --temp_dir TEMP_DIR \             # optional
      --download_species \              # optional
      --dropseq_tools DROPSEQ_TOOLS     # required
welcome to the documentation of spacemake
=========================================

**spacemake: pipeline for processing and analysing sequencing based spatial-transcriptomics data**

.. image:: graphical_abstract_twitter.png
    :width: 800
    :alt: graphical abstract

.. note::

    This project is under active development.

Contents
--------

.. toctree::
    :maxdepth: 3

    install
    quick-start/index.rst
    initialize
    config
    projects/index
    run
    tutorials/index
    api/index

snakemake functions
===================


Installation
============

Step 1: create conda environment
--------------------------------

The most straightforward way of installing spacemake, is first creating a conda environment with the above packages.

To do this, we highly recommend using `mamba <https://github.com/mamba-org/mamba>`_, a much faster conda package manager than conda itself.

After mamba is installed, download the `environment.yaml <https://raw.githubusercontent.com/rajewsky-lab/spacemake/dev/environment.yaml>`_. This environment.yaml, contains all the dependencies required by spacemake.

Once downloaded, to install all spacemake dependencies type::

    mamba env create -f environment.yaml

This will create a conda environment called ``spacemake``.

Too activate the newly created environment type::

   conda activate spacemake

Step 2: download Dropseq-tools
------------------------------

To work with spacemake, currently it is needed to download `Dropseq-tools from here <https://github.com/broadinstitute/Drop-seq>`_.
This packages is a collection of processing tools originally written for `Drop-seq <https://www.cell.com/cell/fulltext/S0092-8674(15)00549-8>`_. Spacemake uses several functions from this package during pre-processing and processing and without it it is impossible to run spacemake.

Simply download one of the releases (we recommend using `2.5.1 <https://github.com/broadinstitute/Drop-seq/releases/download/v2.5.1/Drop-seq_tools-2.5.1.zip>`_) and place it somewhere in your filesystem.


Step 3: install spacemake
-------------------------

**After creating the conda environment and downloading Dropseq-tools** (as described above) spacemake can be installed via ``pip``::

   pip install spacemake

This will install spacemake, you should be good to go :)

.. warning::
    Make sure to first create the conda environment as described above.

    Although it is also possible to install the required packages independently, and then
    to install spacemake, this option has not been tested, and one can quickly run into
    dependency issues and errors.

To install the developmental version of spacemake (``dev`` branch from github) type the following command::

   pip install git+https://github.com/rajewsky-lab/spacemake.git@dev

.. _Running spacemake general:

Running spacemake
=================

Main modules
------------

After spacemake in configured with the ``spacemake config`` command, and projects/samples
are added with the ``spacemake projects`` command, spacemake can be run with the 
``spacemake run`` command. It takes the following parameters::

    spacemake run \ 
        --cores CORES \     # number of cores to be used in total
        --dryrun, -n  \     # invokes a dry snakemake run, printing only commands
        --rerun-incomplete, --ri \
                            # forces snakemake to rerun incompletely generated files
        --keep-going  \     # if a job fails, keep executing independent jobs.
                            # we recommend to always set this when running spacemake
                            # overnight
        --printshellcmds, -p \
                            # print shell commands for each rule, if exist
        --touch, -t   \     # rather than running the rules, just touch each file
        --with_fastqc, -wfqc
                            # Run also fastqc as part of the spacemake run

Downsampling
------------

To run a downsampling (or saturation) analysis, one can use the following command::

    spacemake run downsample \
        --project_id_list [PROJECT_ID_LIST ...] \
        --sample_id_list [SAMPLE_ID_LIST ...]

In the ``project_id_list`` and ``sample_id_list`` arguments one can specify which a 
list of ``project_id``-s and ``sample_id``-s respectively, for which the downsampling
should be run. It is possible to set only one, or both of these arguments. If both are
set the downsampling will be run on samples for which the ``project_id`` and the ``sample_id`` are in both lists (intersection).

.. note::

    In addition to the list arguments specified above, the downsample command also
    takes the same arguments as the simple ``spacemake run`` command.

.. _Seq-scope: https://www.sciencedirect.com/science/article/pii/S0092867421006279
.. _Visium: https://www.10xgenomics.com/products/spatial-gene-expression
.. _Slide-seq: https://www.nature.com/articles/s41587-020-0739-1
.. _Drop-seq: https://mccarrolllab.org/dropseq/
.. _10X Chromium: https://www.10xgenomics.com/products/single-cell-gene-expression
Configuration
=============

Once installed, spacemake configured before running.

.. include:: shared/spacemake_init.rst

Optionally, you can also provide the ``--download_species`` flag, which will download Gencode genomes and
annotations for ``mouse`` and ``human``, and place them under ``project\_root/species\_data/<species>``,
where <species> is either mouse or human.

.. include:: shared/shared_sample_variables.rst

Configure species
-----------------

.. _configure-species:

To add species, the following command can be used::

   spacemake config add_species \
       --name NAME \         # name of the species to be added
       --genome GENOME \     # path to the genome (.fa) file for the species to
                             # be added
       --annotation ANNOTATION \
                             # path to the annotation (.gtf) file for the species
                             # to be added
       --rRNA_genome RRNA_GENOME
                             # (optional) path to the ribosomal-RNA genome (.fa)
                             # file for the species to be added

The ``spacemake config update_species`` takes the same arguments as above, while ``spacemake config delete_species`` takes only ``--name``.

To list the currently available ``species``, type::
   
   spacemake config list_species

Configure barcode\_flavors
--------------------------

.. _configure-barcode_flavor:

This sample-variable describes how the cell-barcode and the UMI should be extracted from Read1 and Read2.
The ``default`` value for barcode\_flavor will be dropseq: ``cell_barcode = r1[0:12]`` (cell-barcode comes from first 12nt of Read1) and
``UMI = r1[12:20]`` (UMI comes from the 13-20 nt of Read1). 

**If a sample has no barcode\_flavor provided, the default run\_mode will be used**

Provided barcode\_flavors
^^^^^^^^^^^^^^^^^^^^^^^^^

Spacemake provides the following barcode\_flavors out of the box:

.. code-block:: yaml

    default:
        cell: "r1[0:12]"
        UMI: "r1[12:20]"
    slide_seq_14bc:
        cell: "r1[0:14]"
        UMI: "r1[14:23]"
    slide_seq_15bc:
        cell: "r1[0:14]"
        UMI: "r1[15:23]"
    visium:
        cell: "r1[0:16]"
        UMI: "r1[16:28]"
    sc_10x_v2:
        cell: "r1[0:16]"
        UMI: "r1[16:26]"
    seq_scope:
        UMI: "r2[0:9]"
        cell: "r1[0:20]"

To list the currently available ``barcode_flavor``-s, type::
   
   spacemake config list_barcode_flavors

Add a new barcode\_flavor
^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block::

   spacemake config add_barcode_flavor \
      --name NAME \
         # name of the barcode flavor

      --umi UMI \
         # structure of UMI, using python's list syntax.
         # Example: to set UMI to 13-20 NT of Read1, use --umi r1[12:20].
         # It is also possible to use the first 8nt of Read2 as UMI: --umi r2[0:8].

      --cell_barcode CELL_BARCODE
         # structure of CELL BARCODE, using python's list syntax.
         # Example: to set the cell_barcode to 1-12 nt of Read1, use --cell_barcode r1[0:12].
         # It is also possible to reverse the CELL BARCODE, for instance with r1[0:12][::-1]. 


Update/delete a barcode\_flavor
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The ``spacemake config update_barcode_flavor`` takes the same arguments as above, while ``spacemake config delete_barcode_flavor`` takes only ``--name``.

Configure run\_modes
--------------------

.. _configure-run_mode:

Specifying a "run mode" is an essential flexibity that spacemake offers.
Through setting a ``run_mode``, a sample can be processed and analysed downstream in various fashions.

Each ``run_mode`` can have the following variables:

``n_beads``
   number of cell-barcode expected

``umi_cutoff``
   a list of integers. downstream the analysis will be run using these UMI cutoffs,
   that is cell-barcodes with less UMIs will be discarded

``clean_dge``
   whether to clean cell-barcodes from overhang primers, before creating the DGE.

``detect_tissue`` (spatial only)
   if ``True``, apart from UMI cutoff spacemake will try to detect the tissue *in-silico*.

``polyA_adapter_trimming``
   if ``True`` 3' polyA stretches and apaters will be trimmed from Read2.

``count_intronic_reads``
   if ``True`` intronic reads will be counted when creating the DGE.

``count_mm_reads``
   if ``True`` multi-mappers will be counted. Only those multi-mapping reads will be
   counted this way, which map to exactly one CDS or UTR segment of a gene.

``mesh_data`` (spatial only)
   if ``True`` a mesh will be created when running this ``run_mode``.

``mesh_type`` (spatial only)
   spacemake currently offers two types of meshes: (1) ``circle``, where circles with a given
   ``mesh_spot_diameter_um`` will be placed in a hexagonal grid, ``mesh_spot_distance_um``
   distance apart; (2) a hexagonal grid, where equal hexagons with ``mesh_spot_diameter_um``
   sides will be placed in a full mesh grid, such that the whole area is covered.

``mesh_spot_diameter_um`` (spatial only)
   the diameter of the mesh spatial-unit, in microns.

``mesh_spot_distance_um`` (spatial only, only for circle mesh)
   distance between the meshed circles, in microns.

``parent_run_mode``
   Each ``run_mode`` can have a parent, to which it will fall back.
   If a one of the ``run_mode`` variables is missing, the variable of the parent will be used.
   If parent is not provided, the ``default`` ``run_mode`` will be the parent. 

Provided run\_mode(s)
^^^^^^^^^^^^^^^^^^^^^

.. code-block:: yaml

    default:
        n_beads: 100000
        umi_cutoff: [100, 300, 500]
        clean_dge: False
        detect_tissue: False
        polyA_adapter_trimming: True
        count_intronic_reads: True
        count_mm_reads: False
        mesh_data: False
        mesh_type: 'circle'
        mesh_spot_diameter_um: 55
        mesh_spot_distance_um: 100
    visium:
        n_beads: 10000
        umi_cutoff: [1000]
        clean_dge: False
        detect_tissue: True
        polyA_adapter_trimming: False
        count_intronic_reads: False
        count_mm_reads: True
    slide_seq:
        n_beads: 100000
        umi_cutoff: [50]
        clean_dge: False
        detect_tissue: False
    scRNA_seq:
        n_beads: 10000
        umi_cutoff: [500]
        detect_tissue: False
        polyA_adapter_trimming: True
        count_intronic_reads: True
        count_mm_reads: False
    seq_scope:
        clean_dge: false
        count_intronic_reads: false
        count_mm_reads: false
        detect_tissue: false
        mesh_data: true
        mesh_spot_diameter_um: 10
        mesh_spot_distance_um: 15
        mesh_type: hexagon
        n_beads: 1000
        umi_cutoff:
        - 100
        - 300

.. note::
   If a sample has no ``run_mode`` provided, the ``default`` will be used

.. note:: 
   If a ``run_mode`` variable is not provided, the variable of the default ``run_mode`` will be used

To list the currently available ``run_mode``-s, type::
   
   spacemake config list_run_modes

Add a new run\_mode
^^^^^^^^^^^^^^^^^^^

See the :ref:`variable descriptions <configure-run_mode>` above.

.. code-block::

   spacemake config add_run_mode \
      --name NAME \ 
      --parent_run_mode PARENT_RUN_MODE \
      --umi_cutoff UMI_CUTOFF [UMI_CUTOFF ...] \
      --n_beads N_BEADS \
      --clean_dge {True,true,False,false} \
      --detect_tissue {True,true,False,false} \
      --polyA_adapter_trimming {True,true,False,false} \
      --count_intronic_reads {True,true,False,false} \
      --count_mm_reads {True,true,False,false} \
      --mesh_data {True,true,False,false} \
      --mesh_type {circle,hexagon} \
      --mesh_spot_diameter_um MESH_SPOT_DIAMETER_UM \
      --mesh_spot_distance_um MESH_SPOT_DISTANCE_UM

Update/delete a run\_mode
^^^^^^^^^^^^^^^^^^^^^^^^^

The ``spacemake config update_run_mode`` takes the same arguments as above, while ``spacemake config delete_run_mode`` takes only ``--name``.


Configure pucks
---------------

.. _configure-puck:

Each spatial sample, needs to have a ``puck``. The ``puck`` sample-variable will define the 
dimensionality of the underlying spatial structure, which then spacemake will use
during the autmated analysis and plotting. 

Each puck has the following variables:

- ``width_um``: the width of the puck, in microns
- ``spot_diameter_um``: the diameter of bead on this puck, in microns.
- ``barcodes`` (optional): the path to the barcode file, containing the cell\_barcode
  and (x,y) position for each. This is handy, when several pucks have the same barcodes,
  such as for 10x visium.


Provided pucks
^^^^^^^^^^^^^^

.. code-block:: yaml

    default:
        width_um: 3000
        spot_diameter_um: 10
    visium:
        barcodes: 'puck_data/visium_barcode_positions.csv'
        width_um: 6500
        spot_diameter_um: 55
    seq_scope:
        width_um: 1000
        spot_diameter_um: 1
    slide_seq:
        width_um: 3000
        spot_diameter_um: 10

as you can see, the ``visium`` puck comes with a ``barcodes`` variable, which points to
``puck_data/visium_barcode_positions.csv``. Upon initiation, this file will automatically placed 
there by spacemake

To list the currently available ``puck``-s, type::
   
   spacemake config list_pucks


Add a new puck
^^^^^^^^^^^^^^

.. code-block::

   spacemake config add_puck \
      --name NAME \        # name of the puck
      --width_um WIDTH_UM \
      --spot_diameter_um SPOT_DIAMETER_UM \
      --barcodes BARCODES # path to the barcode file, optional 

Manage projects and samples
===========================

In spacemake each sample, and it's settings, are stored in the ``project_df.csv`` under the root
directory of the spacemake project.

Each sample will have exactly one row in this ``project_df.csv`` file. In the back-end, spacemake uses a ``pandas.DataFrame`` to load, and save this ``.csv`` file on disk. This data-frame
will be indexed by key ``(project_id, sample_id)``

The spacemake class responsible for this back-end logic is the :ref:`ProjectDF<ProjectDF>` class.

Add a single sample
-------------------

Sample parameters
^^^^^^^^^^^^^^^^^

In spacemake each sample can have the folloing variables:

``project_id``
   ``project_id`` of a sample

``sample_id``
   ``sample_id`` of a sample

``R1``
   ``.fastq.gz`` file path(s) to Read1 read file(s). Can be either a single file, or a space separated list of consecutive files. If a list provided, the files will be merged together and the merged ``R1.fastq.gz`` will be processed downstream.

``R2``
    same as before, but for Read2 read file(s).    

``longreads`` (optional)
   fastq(.gz)|fq(.gz)|bam file path to pacbio long reads for library debugging

``longread-signature`` (optional)
   identify the expected longread signature (see longread.yaml)

``dge`` (optional)
    Since the ``0.1`` version of spacemake, it is possible to only provide the count matrix as input data for spacemake.
    Note: a raw count matrix is expected, if a non count matrix is provided, spacemake will raise an error. 

``barcode_flavor`` (optional)
   ``barcode_flavor`` of the sample. If not provided, ``default`` will be used (Drop-seq).

``species``
   ``species`` of the sample

``puck`` (optional)
   name of the ``puck`` for this sample. if puck contains a ``barcodes`` variable, with a path
   to a coordinate file, those coordinates will be used when processing this sample.
   If not provided, a ``default`` puck will be used with ``width_um=3000``,
   ``spot_diameter_um=10``.

``puck_id`` (optional)
   ``puck_id`` of a sample

``puck_barcode_file`` (optional)
    the path to the file contining (x,y) positions of the barcodes. If the ``puck`` for this
    sample has a ``barcodes`` variable, it will be ignored, and ``puck_barcode_file`` will
    be used.

``investigator`` (optional)
   name investigator(s) responsible for this sample

``experiment`` (optional)
   description of the experiment

``sequencing_date`` (optional)
   sequencing date of the sample

``run_mode`` (optional)
   A list of ``run_mode`` names for this sample. The sample will be processed as defined in 
   the ``run_mode``-s provided. If not provided, the ``default`` ``run_mode`` will be used.


To add a single sample, we can use the following command::

   spacemake projects add_sample \
      --project_id PROJECT_ID \                 # required
      --sample_id SAMPLE_ID \                   # required
      --R1 R1 [R1 R1 ...] \                     # required, if no longreads
      --R2 R2 [R2 R2 ...] \                     # required, if no longreads
      --longreads LONGREADS \                   # required, if no R1 & R2
      --longread-signature LONGREAD_SIGNATURE \ # optional
      --barcode_flavor BARCODE_FLAVOR \         # optional
      --species SPECIES \                       # required
      --puck PUCK \                             # optional
      --puck_id PUCK_ID \                       # optional
      --puck_barcode_file PUCK_BARCODE_FILE \   # optional
      --investigator INVESTIGATOR \             # optional
      --experiment EXPERIMENT \                 # optional
      --sequencing_date SEQUENCING_DATE \       # optional
      --run_mode RUN_MODE [RUN_MODE ...] \      # optional


.. warning::

   A sample is spatial only if: either a ``puck_barcode_file`` is provided, or the sample's
   ``puck`` has a ``barcodes`` variable pointing to a barcode position file.
   If this is not the case, spacemake won't be able to find the spatial barcodes for
   this sample, and the sampe will be processed as a single-cell sample.

   In case both the ``puck_barcode_file`` is provided and the sample's ``puck`` has the
   ``barcodes`` variable set, ``puck_barcode_file`` will be used for the spatial coordinates.

Add a Visium/Seq-scope/Slide-seq sample
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Currently spacemake works out of the box with three spatial methods: `Visium <https://www.10xgenomics.com/products/spatial-gene-expression>`_, `Seq-scope <https://www.sciencedirect.com/science/article/abs/pii/S0092867421006279>`_ and `Slide-seq <https://pubmed.ncbi.nlm.nih.gov/33288904/>`_.

* To add a Visium sample, follow the :ref:`quick start guide here <step 1: add a visium sample>`.
* To add a Seq-scope sample, follow the :ref:`quick start guide here <step 1: add a seq-scope sample>`.
* To add a Slide-seq sample, follow the :ref:`quick start guide here <step 1: add a slide-seq sample>`.

Add a custom spatial sample
^^^^^^^^^^^^^^^^^^^^^^^^^^^

In order to process a custom spatial sample with spacemake follow the step by step guide below.

Step 1: specifying a puck
"""""""""""""""""""""""""

Each spatial sample will need a so-called puck to be configured first. By 'puck' we mean the physical properties of the underlying methods.
Visium for instance works with 6.5mm by 6.5mm sized capture areas, where each spot has 55 microns diameter. To configure a custom puck :ref:`follow the guide here <configure pucks>`.

.. warning::

    If a puck is not specified, spacemake will still run but will use the ``default`` puck as specified :ref:`here <provided pucks>`.

Step 2: formatting a custom puck_barcode_file
"""""""""""""""""""""""""""""""""""""""""""""

For all spatial samples we need to provide a ``puck_barcode_file``. This file needs to be a comma or tab separated, and it needs to have the following three (named) columns:

   - ``cell_bc``, ``barcodes``  or ``barcode`` for cell-barcode
   - ``xcoord`` or ``x_pos`` for x-positions
   - ``ycoord`` or ``y_pos`` for y-positions

Step 3: configure run\_mode(s), barcode\_flavor and species
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Before a custom sample is added the run\_mode(s), barcode\_flavor and species should be configured. The guides on how to do this can be found :ref:`here for run-modes <configure run\\_modes>`, :ref:`here for <configure barcode\\_flavors>` and :ref:`here for species <configure species>`.

The configured run\_mode(s) will specify how a sample is processed downstream, and the barcode\_flavor will specify the barcoding strategy used (ie how many nucleotides are used for UMI, which nucleotides are used for the spot barcodes).

.. warning::

    If no run\_mode(s) are provided spacemake will use the ``default`` run\_mode as specified :ref:`here <provided run\\_mode(s)>`.

    Similarily if there is no barcode\_flavor specified spacemake will use the ``default`` barcode\_flavor as specified :ref:`here <provided barcode\\_flavors>`.

Step 4: add your sample
"""""""""""""""""""""""

Once everything is configured you can add your custom spatial sample with the following command::

    spacemake projects add_sample \
        # your sample's project_id \
        --project_id PROJECT_ID \
        # your sample's sample_id \
        --sample_id SAMPLE_ID \
        # one or more R1.fastq.gz files
        --R1 R1 [R1 R1 ...] \
        # one or more R2.fastq.gz files
        --R2 R2 [R2 R2 ...] \
        # name of the barcode\_flavor, configured in Step 3 \
        --barcode_flavor BARCODE_FLAVOR \
        # name of the species, configured in Step 3 \
        --species SPECIES \
        # name of the puck, configured in Step 1 \
        --puck PUCK \
        # path to your custom barcode file, configured in Step 2 \
        --puck_barcode_file PUCK_BARCODE_FILE \
        # name of the run\_mode(s), configured in Step 3 \
        --run_mode RUN_MODE [RUN_MODE ...]

Add a single-cell sample
^^^^^^^^^^^^^^^^^^^^^^^^

To add a single-cell sample follow the :ref:`quick start guide here <step 1: add a single-cell rna-seq sample>`.

Add a pre-processed count-matrix
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Coming soon!

Add several samples at once
---------------------------

.. _add-several-samples:

It is possible to add several samples in just one command. First, the sample variables have
to be defined in a ``samples.yaml`` file, then we can run the following command::

   spacemake projects add_samples_from_yaml --samples_yaml samples.yaml

The ``samples.yaml`` should have the following structure:

.. code-block:: yaml

   additional_projects:
      - project_id: visium
        sample_id: visium_1
        R1: <path_to_visium_1_R1.fastq.gz>
        R2: <path_to_visium_1_R2.fastq.gz>
        species: mouse
        puck: visium
        barcode_flavor: visium
        run_mode: [visium]
      - project_id: visium
        sample_id: visium_2
        R1: <path_to_visium_2_R1.fastq.gz>
        R2: <path_to_visium_2_R2.fastq.gz>
        species: human
        puck: visium
        barcode_flavor: visium
        run_mode: [visium]
      - project_id: slideseq
        sample_id: slideseq_1
        R1: <path_to_slideseq_1_R1.fastq.gz>
        R2: <path_to_slideseq_1_R2.fastq.gz>
        species: mouse
        puck: slideseq
        barcode_flavor: slideseq_14bc
        run_mode: [default, slideseq]
        puck_barcode_file: <path_to_slideseq_puck_barcode_file>

Under ``additional_projects`` we define a list where each element will be a key:value pair, to be inserted in the ``project_df.csv``

.. note::
   When using the above command, if a sample is already present in the ``project_df.csv`` rather than adding it again, spacemake will update it.
   
   If someone runs ``spacemake projects add_samples_from_yaml --samples yaml samples.yaml`` and
   then modifies something in the ``samples.yaml``, and runs the command again, the ``project_df.csv``
   will contain the updated version of the settings.

Add samples from illumina sample-sheet
--------------------------------------

Coming soon...

Listing projects
----------------

To list projects, which are already configured and added, simply type::
    
    spacemake projects list

It will show the main variables for each project in the ``project_df.csv``. 

To view extra variables which are not shown, use the ``--variables`` option 
to specify which extra variables to show.
One of the most important parts of spacemake are the so-called 'shared sample-variables'.
These are reusable, user-definable variables, which we can assign to several samples.

They can be shortly defined as follows:

species
   a collection of genome, annotation and rRNA\_genome. There is no default species, and each sample can have exactly one species.

barcode\_flavor
   the variable which specifies the structure of Read1 and Read2, namely how the cell\_barcode and UMI should be extracted. If no value provided for a sample, the default will be used.

run\_mode
   each sample can have several ``run_mode``-s, all of which are user definable. If no ``run_mode``-s are specified, a sample will be processed using ``default`` ``run_mode`` settings.

puck (for spatial samples only)
   if a sample is spatial, it has to have a puck variable. If no puck is specified, a default puck will be used.  


To add, update, delete or list a shared sample-variable, you can use the following commands::

   spacemake config add_<shared-sample-variable>
   spacemake config update_<shared-sample-variable>
   spacemake config delete_<shared-sample-variable>
   spacemake config list_<shared-sample-variable>

where ``<shared-sample-variable>`` is one of ``species, barcode_flavor, run_mode or puck``
After you have installed spacemake as specified :ref:`here <installation>`, you are ready to process and analyze spatial samples.

To initialize spacemake ``cd`` into the directory in which you want to start spacemake. This directory will be your ``project_root``.

Then simply type::
   
   spacemake init \
      --dropseq_tools <path_to_dropseq_tools_dir>

Here the `path_to_dropseq_tools_dir` should point to the directory of the downloaded Dropseq-tools package, downloaded :ref:`in Step 2 of the installation <step 2: download dropseq-tools>`.
After a sample is added spacemake can be run with::

   spacemake run --cores <n_cores> --keep-going

The ``--keep-going`` flag is optional, however it will ensure that spacemake runs all
the jobs it can, even if one job fails (this logic is directly taken from snakemake).

For a complete explanation on the `spacemake run` command :ref:`check out the documentation here <Running spacemake general>`.
.. include:: ../links.rst

.. _Quick start guide main:

Quick start guide
=================

The examples here are minimal code pieces on how to quickly get started with spacemake.
It is assumed that spacemake has been instaled following the instructions :ref:`here <installation>`.

.. _Quick start guide initialize spacemake:

Initialize spacemake
--------------------

.. include:: ../shared/spacemake_init.rst


.. _Quick start guide shared sample-variables:

Shared sample-variables
-----------------------

.. include:: ../shared/shared_sample_variables.rst

As spacemake comes with no ``default`` value for ``species``, before anything can be done,
a new species has to be added::

   spacemake config add_species \
      --name \         # name of the species
      --genome \       # path to .fa file
      --annotation \   # path to .gtf file
      --rRNA_genome \  # (optional) path to ribosomal-RNA genome
      --STAR_index_dir # (optional) path to an existing STAR index directory
      
More info :ref:`here <configure-species>`.

.. warning::

    If the ``--STAR_index_dir`` flag is provided spacemake will check if the STAR 
    index provided has the same version of STAR as the command-line STAR. If this is
    not the case, an error will be raised.

Visium quick start
------------------

Step 1: add a Visium sample
^^^^^^^^^^^^^^^^^^^^^^^^^^^

After :ref:`spacemake has been initialized <initialization>`, a `Visium`_ sample can be added.

To add a `Visium`_ sample, type in terminal:

.. code-block:: console

   spacemake projects add_sample \
      --project_id <project_id> \
      --sample_id <sample_id> \
      --R1 <path_to_R1.fastq.gz> \ # single R1 or several R1 files
      --R2 <path_to_R2.fastq.gz> \ # single R2 or several R2 files
      --species <species> \
      --puck visium \
      --run_mode visium \
      --barcode_flavor visium

Above we add a new visium project with ``puck, run_mode, barcode_flavor`` all set to ``visium``.

This is possible as spacemake comes with pre-defined variables, all suited for visium. The visium ``run_mode`` will process the 
sample in the same way as `spaceranger <https://support.10xgenomics.com/spatial-gene-expression/software/pipelines/latest/what-is-space-ranger>`_ would: intronic reads will not be counted, multi-mappers (where the multi-mapping read maps only to one CDS or UTR region) will be counted,
3' polyA stretches will not be trimmed from Read2.

.. note::

   With the ``--R1`` and ``--R2`` it is possible to provide a single ``.fastq.gz`` file (one per mate) or several files per mate.
   For example, if the result of a demultiplexing run is as follows:
   
   ``sample_1a_R1.fastq.gz``, ``sample_1b_R1.fastq.gz``, ``sample_1a_R2.fastq.gz``, ``sample_1b_R2.fastq.gz``, meaning that
   R1 and R2 are both split into two, one can simply call spacemake with the following command::
      
      spacemake projects add_sample \
         ...
         --R1 sample_1a_R1.fastq.gz sample_1b_R1.fastq.gz \
         --R2 sample_1a_R2.fastq.gz sample_1b_R2.fastq.gz \

   The important thing is to always keep the order consistent between the two mates.

To see the values of these predefined variables checkout the :ref:`configuration <Configuration>` docs.

**To add several visium samples at once, follow** :ref:`the tutorial here <add-several-samples>`

.. _running spacemake Visium:

Step 2: running spacemake
^^^^^^^^^^^^^^^^^^^^^^^^^

.. include:: run_spacemake.rst 

Slide-seq quick start
---------------------

Step 1: add a Slide-seq sample
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

After :ref:`spacemake has been initialized <initialization>`, a `Slide-seq`_ sample can be added.

To add a `Slide-seq`_ sample, type in terminal:

.. code-block:: console

   spacemake projects add_sample \
      --project_id <project_id> \
      --sample_id <sample_id> \
      --R1 <path_to_R1.fastq.gz> \
      --R2 <path_to_R2.fastq.gz> \
      --species <species> \
      --puck slide_seq \
      --run_mode slide_seq \
      --barcode_flavor slide_seq_14bc \
      --puck_barcode_file <path_to_puck_barcode_file>

Above we add a new `Slide-seq`_ project with the ``puck, run_mode`` will be set to ``slide_seq``
which are pre-defined settings for `Slide-seq`_ samples.

.. note::
   For spatial samples other than visium - such as `Slide-seq`_ - we need to provide a
   ``puck_barcode_file`` (since each puck has different barcodes, unlike for visium samples).
   This file should be a comma or tab separated, containing column names as first row. Acceptable column names are:

   - ``cell_bc``, ``barcodes``  or ``barcode`` for cell-barcode
   - ``xcoord`` or ``x_pos`` for x-positions
   - ``ycoord`` or ``y_pos`` for y-positions

In this example ``barcode_flavor`` will be set to ``slide_seq_14bc``,
a pre-defined ``barcode_flavor`` in spacemake, where the ``cell_barcode`` comes from the first 14nt of Read1, and the ``UMI`` comes from nt 13-22 (remaining 9 nt). 
The other pre-defined ``barcode_flavor`` for `Slide-seq`_ is ``slide_seq_15bc``: here ``cell_barcode`` again comes from the first 14nt of Read1, but the ``UMI`` comes from nt 14-22 (remaining 8) of Read1.

To see the values of these predefined variables checkout the :ref:`configuration <Configuration>` docs.

**To add several slide_seq projects at once, follow** :ref:`the tutorial here <add-several-samples>`

.. _running spacemake Slide-seq:

Step 2: running spacemake
^^^^^^^^^^^^^^^^^^^^^^^^^

.. include:: run_spacemake.rst 

Seq-scope quick start
---------------------

Step 1: add a Seq-scope sample
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

After :ref:`spacemake has been initialized <initialization>`, a `Seq-scope`_ sample can be added.

Adding a `Seq-scope`_ sample 
is similar to Slide-seq:

.. code-block:: console

   spacemake projects add_sample \
      --project_id <project_id> \
      --sample_id <sample_id> \
      --R1 <path_to_R1.fastq.gz> \ # single R1 or several R1 files
      --R2 <path_to_R2.fastq.gz> \ # single R2 or several R2 files
      --species <species> \
      --puck seq_scope \
      --run_mode seq_scope \
      --barcode_flavor seq_scope \
      --puck_barcode_file <path_to_puck_barcode_file>

Here we used the pre-defined variables for ``puck``, ``barcode_flavor`` and ``run_mode`` all set to ``seq_scope``.

The ``seq_scope`` ``puck`` has 1000 micron width and bead size set to 1 micron.

The ``seq_scope`` ``barcode_flavor`` describes how the ``cell_barcode`` and he ``UMI`` should be extracted. 
As described in the `Seq-scope`_ paper.
``cell_barcode`` comes from nt 1-20 of Read1, and ``UMI`` comes from 1-9nt of Read2.

The ``seq_scope`` ``run_mode`` has its settings as follows:

.. code-block:: yaml

    seq_scope:
        clean_dge: false
        count_intronic_reads: false
        count_mm_reads: false
        detect_tissue: false
        mesh_data: true
        mesh_spot_diameter_um: 10
        mesh_type: hexagon
        n_beads: 1000
        umi_cutoff:
        - 100
        - 300

The most important thing to notice here that by default, we create a hexagonal mesh with
the ``seq_scope`` ``run_mode``. This means that downstream rather than with working with
the 1 micron beads, spaceame will create a mesh of adjascent, equal hexagons with 10 micron
sides.

.. _running spacemake Seq-scope:

Step 2: running spacemake
^^^^^^^^^^^^^^^^^^^^^^^^^

.. include:: run_spacemake.rst 

scRNA-seq quick start
---------------------

Step 1: add a single-cell RNA-seq sample
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

spacemake was written as a spatial-transcriptomics pipeline, however it will also work for
single-cell experiments, where there is no spatial information available. 

To add a scRNA-seq sample, simply type::

   spacemake projects add_sample \
      --project_id <project_id> \
      --sample_id <sample_id> \
      --R1 <path_to_R1.fastq.gz> \ # single R1 or several R1 files
      --R2 <path_to_R2.fastq.gz> \ # single R2 or several R2 files
      --species <species> \
      --run_mode scRNA_seq \
      --barcode_flavor default # use other flavors for 10X Chromium

As seen above, we define fewer variables as before: only ``species`` and ``run_mode`` are needed.

.. warning:: 

    As it can be seen, we set the ``barcode_flavor`` to ``default``, which will use the `Drop-seq`_ 
    scRNA-seq barcoding strategy (``cell_barcode`` is 1-12nt of Read1, ``UMI`` is 13-20 nt of Read1).

    For 10X samples either the ``sc_10x_v2`` (10X Chromium Single Cell 3' V2)
    or ``visium`` (10X Chromium Single Cell 3' V3, same as visium) should be be used as
    ``barcode_flavor``. Both are pre-defined in spacemake.
    
    If another barcoding strategy is used, a custom ``barcode_flavor`` has to be defined.
    :ref:`Here is a complete guide on how to add a custom barcode_flavor <add a new barcode\\_flavor>`.

.. note:: 

    By setting ``run_mode`` to ``scRNA_seq`` we used the pre-defined ``run_mode`` settings tailored for single-cell experiments: expected number of cells (or beads) will be 10k, introns will be counted, UMI cutoff will be at 500, multi-mappers will not be counted and polyA and adapter sequences will be trimmed from Read2. 

    Of course, running single-cell samples with other ``run_mode`` settings is also possible.

    :ref:`Here it can be learned how to add a custom run_mode <add a new run\\_mode>`, tailored to the experimental needs.

To see the values of these predefined variables checkout the :ref:`configuration <Configuration>` docs.

**To add several single-cell projects at once, follow** :ref:`the tutorial here <add-several-samples>`

.. _running spacemake scRNA-seq:

Step 2: running spacemake
^^^^^^^^^^^^^^^^^^^^^^^^^

.. include:: run_spacemake.rst 

Longreads Integration
=====================

A Spacemake sample can not only be associated with Illumina sequencing reads 
via the ``--R1`` and ``--R2`` command line arguments, but can also be assigned 
long reads (e.g. PacBio, Nanopore, etc.) using the ``--longreads`` command-line 
argument.

We have added this functionality for the purpose of trouble-shooting 
problems with library construction (see the Spacemake paper). Long sequencing reads can
capture the entire cDNA sequence, including Illumina sequencing adapters, primer handles, 
cell and UMI barcodes, polyA/T as well as the actual cDNA insert.

If you add long reads to a sample, you will enable spacemake to annotate a catalog of 
oligo sequences, thought of as building blocks of your library, against every long read. 
This allows to assess the following feaures:

  - the fraction of molecules which conform to the expected layout of building blocks.
    For example, SMART-handle, (barcode), polyT, (cDNA insert), TSO-handle. Parts in 
    parenthesis here are inferred from the spacing between the adjacent blocks but not 
    directly annotated

  - the fraction of molecules which are missing one or more the expected building blocks
  
  - size distributions of all matches and distributions of their start and end positions,
    which in turn allows to infer the size distributions of inserts and barcode sequences.

  - *concatamerizations* and unexpected, *multiple* occurrences of any building block, 
  - pointing to undesired side-reactions.


Building blocks and Signatures
------------------------------

In order to fully utilize the longread functionality, spacemake needs to know about the sequences 
of all building blocks, as well as their expected layout. 

Tying it together is a ``signature``. Signatures are defined (at the moment manually) in 
a ``longread.yaml`` file. If spacemake does not find a file with that name at the root of 
your ``projects`` folder, it will default to loading ``<spacemake-root>/spacemake/data/longread.yaml``.

Here is an excerpt from this file listing known sequence building blocks.

.. code-block:: yaml

    blocks:
        P5: AATGATACGGCGACCACCGAGATCTACACGCCTGTCCGCGG
        N70X: CTGTCTCTTATACACATCTCCGAGCCCACGAGACNNNNNNNNATCTCGTATGCCGTCTTCTGCTTG
        SMART_primer: AAGCAGTGGTATCAACGCAGAGT
        SMART_bead: AAGCAGTGGTATCAACGCAGAGTAC
        dN-SMRT: AAGCAGTGGTATCAACGCAGAGTGA
        TSO: AAGCAGTGGTATCAACGCAGAGTGAATGGG
        polyT: TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT


You can see here that ``SMART_primer`` and ``SMART_bead`` are very similar sequences, 
differing only at the end, with ``SMART_bead`` containing an extra ``AC`` dinucleotide. 
In case of multiple (partial) matches to a building block, a match with higher score 
(ie more matches) will supersede an overlapping match with lower score. This allows to
distinguish the handle at the start of the capture oligos (which features the ``AC`` at the 3' end) 
from other occurrences of the SMART primer site, which might be introduced via a different route.

The next section in the YAML file describes the expected structure of sequencing reads in terms of
the building blocks. Here is an excerpt including the definition for dropseq/Nadia cDNA libraries:

.. code-block:: yaml

    signatures:
        chromium:
            label: chromium
            prio: 1
            color: gray
            CB: r1[0:16]
            UMI: r1[16:26]
            intact: 10X_start,polyT,10X_TSO_RC
            other: 10X_C3_RT_PRIMER,10X_C2_RT_PRIMER
            cDNA_after: polyT
            prefixes: P5
            suffixes: N70X
            read1_primer: 10X_start
            read2_primer: 10X_TSO

        dropseq:
            label: dropseq
            prio: 2
            color: gray
            CB: r1[8:20]
            UMI: r1[0:8]
            intact: SMART_bead,polyT
            cDNA_after: polyT
            other: SMART_primer,dN-SMRT,TSO,sc_primer
            prefixes: P5
            suffixes: N70X
            read1_primer: SMART_bead
            read2_primer: N70X


Other pre-defined signatures at the moment are ``visium`` (almost identical to chromium) and a couple of
experimental designs from ongoing in-house developments.

The ``prio`` and ``color`` fields are only relevant for overview plots across multiple samples 
and will affect the ordering (prio) and color for the visual representation of sample metrics.

Most important is the ``intact`` field, which lists all buidling blocks expected to be present (in that order) 
on a correct library molecule. P5 and N70X are considered optional prefixes/suffixes at this point, 
because you may choose to perform long read sequencing on a library before or 
after index PCR.

The first building block listed in ``intact`` is expected to be present in all
molecules that derive from the used capture technology (and not some contamination or artifact).
In the case of dropseq/nadia beads, this would be the SMART handle, followed by ``AC``. 
In the case of 10X Chromium or Visium libraries, it would be the ``10X_start`` primer handle attached to
the gel beads or visium slide, respectively.

Occurrences of this first building block are used to distinguish captured molecules from 'other' and 
to orient every long read (sequenced long reads can be either in forward, or reverse-complement orientation).

What happens?
-------------

As soon as at least one of your samples is associated with long sequencing reads, ``spacemake run`` 
will invoke some dedicated tools. Specifically

   1. long reads will be *aligned* against all known building blocks
   2. (overlapping) matches will be *collected* and integrated for each read
   3. based on the presence/absence of each block, each read will be *classified*
   4. *statistics* on the observed long read classes will be gathered, with particular emphasis on the 
      reads falling into the class defined as ``intact`` .
   5. *cDNA* will be extracted and mapped to the genome via `STAR-long`
   6. association of mappability and building block presence/absence is investigated
   7. report *plots* are generated for each sample in ``/processed_data/{sample_id}/longread/reports``

After these steps are completed for every sample with long reads, *overview plots* are generated, which 
present high level results across all samples, side-by-side in ``<projects folder>/longread_overview``.

If you like to utilize any of these functions outside of the spacemake/snakemake workflow you can either 
invoke the longread command via ``python -m spacemake.longread`` or by importing the ``spacemake.longread``
module from your own python scripts.

A nanopre example
-----------------

Here is a full example using a small test data-set. We will download the test data, add the 
sample to a spacemake project, run the analysis, and have a look at the output generated. We use ``wget`` 
to download a small test data set. Alternatively, you can use 
``git clone https://github.com/rajewsky-lab/spacemake-test-data.git`` to check out a collection of different
test data. Here, we further assume that a species labeled ``human`` has already been set up for spacemake [TODO: add link].

.. code-block:: console

   wget https://bimsbstatic.mdc-berlin.de/rajewsky/spacemake-test-data/longread/SRR9008425_subsample.fastq.gz

   spacemake projects add_sample --project_id=test --sample_id=test_longread --longreads=SRR9008425_subsample.fastq.gz --longread-signature=chromium --species=human
   spacemake run -np --cores=64

Spacemake lays out the following tasks to process our longread sample (shortened for brevity):

.. code-block:: console

    Job counts:
            count   jobs
            1       all
            1       cmd_align
            1       cmd_alnstats
            1       cmd_annotate
            1       cmd_edits
            1       cmd_extract
            1       cmd_overview
            1       cmd_report
            1       map_cDNA
            9
    This was a dry-run (flag -n). The order of jobs does not reflect the order of execution.

    rule cmd_align:
        python -m spacemake.longread --parallel=64 --config=longread.yaml --cache=projects/test/processed_data/test_longread/longread/cache/ --annotation-out=projects/test/processed_data/test_longread/longread/annotation/ --stats-out=projects/test/processed_data/test_longread/longread/stats/ --report-out=projects/test/processed_data/test_longread/longread/reports/ --examples-out=projects/test/processed_data/test_longread/longread/examples/ --sample=test_longread --signature=chromium align SRR9008425_subsample.fastq.gz 

    rule cmd_extract:
        python -m spacemake.longread --parallel=1 --config=longread.yaml --cache=projects/test/processed_data/test_longread/longread/cache/ --annotation-out=projects/test/processed_data/test_longread/longread/annotation/ --stats-out=projects/test/processed_data/test_longread/longread/stats/ --report-out=projects/test/processed_data/test_longread/longread/reports/ --examples-out=projects/test/processed_data/test_longread/longread/examples/ --sample=test_longread --signature=chromium extract SRR9008425_subsample.fastq.gz 2> projects/test/processed_data/test_longread/longread/cDNA/test_longread.log > projects/test/processed_data/test_longread/longread/cDNA/test_longread.fa

    rule cmd_annotate:
        python -m spacemake.longread --parallel=1 --config=longread.yaml --cache=projects/test/processed_data/test_longread/longread/cache/ --annotation-out=projects/test/processed_data/test_longread/longread/annotation/     --stats-out=projects/test/processed_data/test_longread/longread/stats/ --report-out=projects/test/processed_data/test_longread/longread/reports/     --examples-out=projects/test/processed_data/test_longread/longread/examples/ --sample=test_longread --signature=chromium annotate SRR9008425_subsample.fastq.gz
    
    rule map_cDNA:
        mkdir -p projects/test/processed_data/test_longread/longread/cDNA/tmp/
        STARlong --runThreadN 8 --genomeDir species_data/human/star_index --genomeLoad NoSharedMemory --readFilesIn projects/test/processed_data/test_longread/longread/cDNA/test_longread.fa --readFilesType Fastx --outSAMtype BAM Unsorted --outSAMunmapped Within --outSAMattributes All --outSAMprimaryFlag AllBestScore --outStd BAM_Unsorted --outFilterMultimapScoreRange 2 --outFilterScoreMin 0 --outFilterScoreMinOverLread 0 --outFilterMatchNminOverLread 0 --outFilterMatchNmin 30 --outFilterMismatchNmax 1000 --winAnchorMultimapNmax 200 --seedSearchStartLmax 12 --seedPerReadNmax 100000 --seedPerWindowNmax 100 --alignTranscriptsPerReadNmax 100000 --alignTranscriptsPerWindowNmax 10000 --outFileNamePrefix projects/test/processed_data/test_longread/longread/cDNA/tmp | /data/rajewsky/shared_bins/Drop-seq_tools-2.4.0/TagReadWithGeneFunction I=/dev/stdin O=projects/test/processed_data/test_longread/longread/cDNA/test_longread.bam ANNOTATIONS_FILE=species_data/human/annotation.gtf

    rule cmd_report:
        python -m spacemake.longread --parallel=1 --config=longread.yaml --cache=projects/test/processed_data/test_longread/longread/cache/     --annotation-out=projects/test/processed_data/test_longread/longread/annotation/     --stats-out=projects/test/processed_data/test_longread/longread/stats/     --report-out=projects/test/processed_data/test_longread/longread/reports/     --examples-out=projects/test/processed_data/test_longread/longread/examples/ --sample=test_longread --signature=chromium report

    rule cmd_edits:
        python -m spacemake.longread --parallel=1 --config=longread.yaml --cache=projects/test/processed_data/test_longread/longread/cache/     --annotation-out=projects/test/processed_data/test_longread/longread/annotation/     --stats-out=projects/test/processed_data/test_longread/longread/stats/     --report-out=projects/test/processed_data/test_longread/longread/reports/     --examples-out=projects/test/processed_data/test_longread/longread/examples/     --sample=test_longread     --signature=chromium       edits SRR9008425_subsample.fastq.gz

    rule cmd_overview:
        python -m spacemake.longread --parallel=1 --config=longread.yaml overview --output longread_overview/ projects/test/processed_data/test_longread/longread/stats/test_longread.report.tsv


Ok, sounds good. Let's do this.

.. code-block:: console

    spacemake run -p --cores=64

.. note:: 

    please adapt ``--cores=64`` to a number appropriate for your machine. The `align` stage of the longread processing is slow
    and benefits strongly from parallelization, we therefore recommend running on a machine with many cores.


Graphical reports
-----------------

After `spacemake run` has finished, you can check the `longread/reports` directory under each sample that was assigned some longread data. You will find several graphical reports
here. Let's walk through the example nanopore data:

Library overview plot
^^^^^^^^^^^^^^^^^^^^^

.. image:: img/test_longread.donuts.png
    :width: 100%

The three panels in the `test_longread.donuts.pdf` plot provide an overview of the library. The donut on the top-left
shows that the majority of reads contain at least ``10X_start`` oligo matches, which means they derive from the capture technology and are labeled `bead-related`.
This, of course is a good thing. The remaining slices of the pie-chart correspond to reads where the 10X_start sequence was either absent or mutated enough to match more with other, related building-blocks. 
A certain fraction of such reads is to be expected and normal. 

The donut plot on the right is a zoom into the `bead-related` section of the left donut plot. Here, we can see that slightly less than half of all reads conform to the expected, full signature. About 50% lack a recognizable match to the template-switch oligo, which could either indicate incomplete long reads or a nanopore library construction strategy which truncates the TSO.

Another useful way to look at the data is the horizontal bar-plot on the lower left. Please note that the x-axis here is logarithmic. Consistent with what we just saw, the most abundant species of long reads is classified as `10X_start,polyT`, i.e. lacking the TSO match. However, this is closely followed by `10X_start,polyT,10X_TSO_RC`, which is the expected species. Note that the `_RC` at the end indicates the reverse complement, which is indeed expected here.
Together these two read species account for the vast majority of the library.

Zooming in on the dependence between the expected building blocks is the last plot, on the lower left. The first bar corresponds to the green `bead-related` slice in the first donut plot. If there were no artifacts or contaminations, one could read this as the fraction of capture oligos on the gel beads which were correctly synthesized such as to begin with `10X_start`.
The next bar-plot shows the fraction of all reads, *starting with `10X_start`* which then go on to display a polyT match next, i.e. they capture the dependence between `10X_start` and `polyT`. As you can see, this is strongly coupled with much more than 95%. In contrast, the last bar is again a condensation of our earlier finding that many reads lack a good TSO match. It represents the fraction of all reads which begin with `10X_start,polyT`, which then continue to contain `10X_TSO_RC`.


Oligo alignment quality plot
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. image:: img/test_longread.oligo_edits.png
    :width: 100%

In order to troubleshoot issues with the building block sequences, up to 1,000 instances of detected building block matches are sampled from the reads and their alignments against the reference sequence is condensed as a step-plot (left side) with reference sequence position on the x-axis and per-base match rate on the y-axis.
To provide more detail, the right side details the kinds of mismatches , deletion and insertions which are observed (y-axis) at each position (x-axis) as a heatmap. Each row in this plot corresponds to one of the expected building blocks listed in ``intact``.

Here, we finally can find a clue as to why TSO-matches are underrepresented: The second half of the TSO sequence is very rarely part of the sequencing reads, perhaps indicating that the nanopore sequencing library was generated with a nested primer, masking part of the TSO, or we are not aligning with the correct reference sequence. However, the heatmap indicates neither mismatches nor insertions, only deletions from position 12 of the TSO_RC sequence on, more consistent with the idea that the reads "just end" before the entire TSO sequence is read.


Histograms
^^^^^^^^^^

.. image:: img/test_longread.hists.png
    :width: 100%

Lastly, the longread module generates cumulative histograms for the start positions (first row), end-positions (second row), match scores (third row) and match lengts (fourth row) of the expected building blocks (columns). Importantly, this plot distinguishes between occurrences that are part of `intact` signature reads and `all` occurrences of a given building block. This allows to quickly spot if we are looking at perhaps different species of molecules/artifacts which are associated with unexpected positioning or spacing of the blocks. This is not the case here, again consistent with the idea that the TSO matches are simply truncated and reads otherwise conforming to expectations.Tutorials
=========

.. toctree::
    he_integration
    manual_he_integration
    process_single_cell_data
    longreads
.. include:: ../links.rst

Processing a custom single-cell sample
======================================

In this tutorial we will process a custom single cell sample. 

As an example we will be using 1 million reads from `this Visium dataset <https://www.10xgenomics.com/resources/datasets/mouse-brain-section-coronal-1-standard-1-0-0>`_.

.. note::
    
    Firstly, the example data used here is a 10X `Visium`_ dataset, hence it is spatial.
    However, for the sake of this tutorial, we will be treating it as a single-cell sample.
    
    Secondly, for many methods (such as `Visium`_, `10X Chromium`_ `Slide-seq`_ or `Seq-scope`_)
    spacemake provides pre-defined variables. If you are using
    one of these methods follow our :ref:`Quick start guide <quick start guide>` instead.

Step 1: install and initialize spacemake
-----------------------------------------

To install spacemake follow the :ref:`installation guide here <installation>`.

To initialize spacemake follow the :ref:`initialization guide here <initialization>`.

Step 2: download test data
--------------------------

For the sake of this tutorial we will work with a test dataset: 1 million Read1 and 1 million Read2 reads from a `Visium`_ adult mouse brain.

To download the test data:

.. code-block::

    wget -nv http://bimsbstatic.mdc-berlin.de/rajewsky/spacemake-test-data/visium/test_fastq/visium_public_lane_joined_1m_R1.fastq.gz
    wget -nv http://bimsbstatic.mdc-berlin.de/rajewsky/spacemake-test-data/visium/test_fastq/visium_public_lane_joined_1m_R2.fastq.gz

.. note:: 

    If there is already data available, to be processed and analyzed, this step can be omitted.

Step 3: add a new species
-------------------------

.. note::

    If you initialized spacemake with the ``--download-species`` flag, you can
    omit this step, as spacemake will automatically download and configure
    mm10 mouse genome.fa and annotation.gtf files for you.

The sample we are working with here is a mouse brain sample, so we have to add a new species:

.. code-block:: console

   spacemake config add_species --name mouse \
   --annotation /path/to/mouse/annotation.gtf \
   --genome /path/to/mouse/genome.fa


Step 4: add a new barcode\_flavor
---------------------------------

The ``barcode_flavor`` will decide which nucletodies of Read1/Read2 extract the UMIs and cell-barcodes from.

In this perticular test sample, the first 16 nucleotides of Read1 are the cell-barcode, and the following 12 nucleotides are the UMIs.

Consequently, we create a new ``barcode_flavor`` like this:

.. code-block:: console

    spacemake config add_barcode_flavor --name test_barcode_flavor \
    --cell_barcode r1[0:16] \
    --umi r1[16:28]

.. note:: 

    There are several ``barcode_flavors`` provided by spacemake out of the box,
    such as ``visium`` for 10X `Visium`_ or ``sc_10x_v2`` for `10X Chromium`_ v2 
    kits. The ``default`` flavor is identical to a `Drop-seq`_ library, with 12
    nucleotide cell-barcode and 8 nucleotide UMI. 

    :ref:`More info about provided flavors here <provided barcode\\_flavors>`.

    If you want to use one of these, there is no need to add your own flavor.

Step 5: add a new run\_mode
---------------------------

A ``run_mode`` in spacemake defines how a sample should processed downstream. 
In this tutorial, we will trim the PolyA stretches from the 3' end of Read2,
count both exonic and intronic reads, expect 5000 cells, and analyze the data,
turn off multi-mapper counting (so only unique reads are counted),
using 50, 100 and 300 UMI cutoffs. To set these parameters, we define a 
``test_run_mode`` like this:

.. code-block:: console

    spacemake config add_run_mode --name test_run_mode \
    --polyA_adapter_trimming True \
    --count_mm_reads False \
    --n_beads 5000 \
    --count_intronic_reads True \
    --umi_cutoff 50 100 300

.. note:: 

    As with ``barcode_flavors``, spacemake provides several ``run_modes`` out
    of the box. For more info :ref:`check out a more detailed guide here <configure run\\_modes>`.

Step 6: add the sample
----------------------

After configuring all the steps above, we are ready to add our (test) sample:

.. code-block:: console

    spacemake projects add_sample --project_id test_project \
    --sample_id test_sample \
    --R1 visium_public_lane_joined_1m_R1.fastq.gz \
    --R2 visium_public_lane_joined_1m_R1.fastq.gz \
    --species mouse \
    --barcode_flavor test_barcode_flavor \
    --run_mode test_run_mode

.. note::

    If there is already data available, here the Read1 and Read2 ``.fastq.gz`` files should be added,
    instead of the test files.

Step 7: runn spacemake
----------------------

Now we can process our samples with spacemake. Since we added only one sample, only one sample will be processed
and analyzed. To start spacemake, simply write:

.. code-block:: console
    
    spacemake run --cores 16

.. note::
    
    The number of cores used should be suited for the machine on which spacemake is ran.
    When processing more than one samle, we recommend using spacemake with at least 8 cores.
    In order to achieve maximum parallelism.

Step 8: results 
---------------

The results of the analysis for this sample will be under ``projects/test_project/processed_data/test_sample/illumina/complete_data/``

Under this directory, there are several files and directories which are important:

* ``final.polyA_adapter_trimmed.bam``: final, mapped, tagged ``.bam`` file. ``CB`` tag contains the cell barcode, and the ``MI`` contains the UMI-s. 

* ``qc_sheet_test_sample_no_spatial_data.html``: the QC-sheet for this sample, as a self-contained ``.html`` file.

* ``dge/``: a directory containing the Digital Expression Matrices (DGEs)

    * ``dge.all.polyA_adapter_trimmed.5000_beads.txt.gz``: a compressed, text based DGE

    * ``dge.all.polyA_adapter_trimmed.5000_beads.h5ad``: the same DGE but stored in ``.h5ad`` format (`used by the anndata python package <https://github.com/theislab/anndata/issues/180>`_). This matrix is stored as a Compressed Sparse Column matrix (using `scipy.sparse.csc_matrix <https://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.csc_matrix.html>`_).

    * ``dge.all.polyA_adapter_trimmed.5000_beads.summary.txt``: the summary of the DGE, one line per cell.

    * ``dge.all.polyA_adapter_trimmed.5000_beads.obs.csv``: the observation table of the matrix. Similar to the previous file, more detailed.

* ``automated_analysis/test_run_mode/umi_cutoff_50/``: In this directory the results of the automated analysis can be found. As it can be seen under the ``automated_analysis`` directory there are two further levels, one for ``run_mode`` and one for ``umi_cutoff``. This is because one sample can have several ``run_modes`` and in the same way one ``run_mode`` can have several UMI cutoffs.

    * ``results.h5ad``: the result of the automated analysis, stored in an anndata object. Same as the DGE before, but containing processed data.

    * ``test_sample_no_spatial_data_illumina_automated_report.html``: automated analysis self-contained ``.html`` report.

.. note::

    If the ``test_project`` had more samples, than those would be automatically placed under ``projects/test_project``. Similarily, under one spacemake
    directory there can be several projects in parallel, and each will have their own directory structure under the ``projects/`` folder.

Manual H&E alignment
====================

Before you begin
----------------

Before start, make sure that you have installed spacemake as specified :ref:`here <installation>`

For the manual allignment we will use Fiji, an open-source image processing tool. Download it from `here <https://imagej.net/software/fiji/downloads>`_.

We will be using tile nr 2105 from `Seq-scope <https://www.sciencedirect.com/science/article/abs/pii/S0092867421006279>`_ for this tutorial. The corresponding H&E image is
`wt_4X_2.jpg <https://deepblue.lib.umich.edu/data/downloads/qv33rw833>`_.

Step 1 - generate an expression image
-------------------------------------

First using the command line we generate an aggregated expression image. In the directory of your spacemake project, type:

.. code-block:: console

    spacemake spatial create_aggregated_expression_img \
        --project_id seq_scope \
        --sample_id seq_scope_liver_2105 \
        --run_mode seq_scope \
        --processed_data False \
        --binary True \
        --out_img aggregated_seq_scope_2105_img_bw.png

This will generate a black and white image based on expression data.

Step 2 - load images into Fiji
------------------------------

In the next step we load both images into Fiji like below:

.. image:: img/manual_alignment_1.png
    :width: 100%
    :alt: Manual alignment first step

Step 3 - select corresponding points
------------------------------------

Next, using the *Multi-point Tool* we manually select corresponding points between our expression image and the H&E image. 
Select a point on one of the images, and then select a corresponding point on the other image. Do this for at least 4-5 corresponding points for a better match.

.. image:: img/manual_alignment_2.png
    :width: 100%
    :alt: Manual alignment second step

Step 4 - align the images
-------------------------

We then use the `Landmark Correspondences <https://imagej.net/plugins/landmark-correspondences>`_ plugin to align the two images based on the correspondencing points we
selected in the previous step. We go to *Plugins -> Transform -> Landmark Correspondences*:

.. image:: img/manual_alignment_3.png
    :width: 100%

In the pop-up window we select H&E image as the *source image* and expression image as the *template image*.
For the *transformation method* select *Moving Least Squares (non-linear)*. Set the *alpha* to *1.00* and the *mesh resolution* to *32*.
Set the *transformation class* to *Affine*.

.. image:: img/manual_alignment_4.png
    :width: 100%

After the transformation we have the two images aligned. We can now save our transformed H&E image (which is aligned with our spatial data).

.. image:: img/manual_alignment_5.png
    :width: 100%


Step 5 - attach the aligned image
---------------------------------

First we load the spacemake processed Seq-scope tile nr 2105 data:

.. code-block:: ipython3

    from spacemake import Spacemake

    spmk = Spacemake('/path/to/your/spacemake/project')

    adata_2105 = spmk.load_processed_adata(
        project_id = 'seq_scope',
        sample_id = 'seq_scope_liver_2105',
        run_mode_name = 'seq_scope',
        umi_cutoff = 300
    )

Then we load the previously manually aligned image and attach it to our data:

.. code-block:: ipython3

    from spacemake.spatial.he_integration import attach_he_adata
    import cv2

    matched_he = cv2.imread('./Transformedwt_4X_2.tif')

    adata = attach_he_adata(adata_2105.copy(),
                            matched_he,
                            push_by_spot_diameter=False,
                            raw_aligned=True)

After attachment, we can plot our expression data on top of the aligned H&E with `scanpy <https://github.com/theislab/scanpy>`_:

.. code-block:: ipython3

    import scanpy as sc

    sc.set_figure_params(dpi=300)

    sc.pl.spatial(adata, color='total_counts')

.. image:: img/manual_alignment_6.png
    :width: 100%


.. note::
    
    As axes in scanpy are flipped with respect to the axes in Fiji, because Fiji reads the image axes in different order.
API
===

Spacemake class
---------------

Accessing spacemake objects from python

.. autoclass:: spacemake.Spacemake
    :members:

H&E integration module
----------------------

.. autofunction:: spacemake.spatial.he_integration.align_he_spot_img

.. autofunction:: spacemake.spatial.he_integration.align_he_aggregated_img

.. autofunction:: spacemake.spatial.he_integration.attach_he_adata

novosparc integration module
----------------------------

.. autofunction:: spacemake.spatial.novosparc_integration.novosparc_denovo

.. autofunction:: spacemake.spatial.novosparc_integration.save_novosparc_res

.. autofunction:: spacemake.spatial.novosparc_integration.novosparc_mapping

.. autofunction:: spacemake.spatial.novosparc_integration.quantify_clusters_spatially
API and Internal API
====================

.. toctree::
    
    api
    internal_api    
Internal API
============

ProjectDF
---------

The ProjectDF class is the core back-end class of spacemake. 

.. autoclass:: spacemake.project_df.ProjectDF
   :members:

ConfigFile
----------

This class is responsible for updating spacemake's configuration.

.. autoclass:: spacemake.config.ConfigFile
   :members:
