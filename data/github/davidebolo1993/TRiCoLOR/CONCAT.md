[![Build Status](https://travis-ci.com/davidebolo1993/TRiCoLOR.svg?branch=master)](https://travis-ci.com/davidebolo1993/TRiCoLOR)

# TRiCoLOR

![alt text](TRiCoLOR.png)

## TRiCoLOR: tandem repeats profiler for long-read sequencing data

TRiCoLOR is available to be installed from source or as a [Docker container](https://hub.docker.com/r/davidebolo1993/tricolor).
Please have a look at [TRiCoLOR documentation](https://davidebolo1993.github.io/tricolordoc/) for any installation or usage questions.

[Source Code](https://github.com/davidebolo1993/TRiCoLOR/tree/master/TRiCoLOR)

[Documentation](https://davidebolo1993.github.io/tricolordoc/)

## Citation

Are you using TRiCoLOR in your works? Please cite:

> Davide Bolognini, Alberto Magi, Vladimir Benes, Jan O Korbel, Tobias Rausch. [TRiCoLOR: tandem repeat profiling using whole-genome long-read sequencing data](https://academic.oup.com/gigascience/article/9/10/giaa101/5918863?searchresult=1). GigaScience. 2020 Oct 7.
These folders contain BAM files simulated with [VISOR](https://academic.oup.com/bioinformatics/article/doi/10.1093/bioinformatics/btz719/5582674/) (accuracy ~ 0.9, reads length ~ 8000 bps, substitutions:insertions:deletions rate ~ 45:25:30, coverage ~ 50X).
BAM files for SON and PARENT1 contain a heterozigous TR expansion.
Normal TR is chr20:17553794-17553824 (15xAC), taken from the list of TRs availables for the GRCh38 human reference genome.
Expanded TR is 115XAC (expansion of 100 motifs).
TRiCoLOR can be applied to identify and profile this expansion in both child and parents as explained in the [manual](https://davidebolo1993.github.io/tricolordoc/).
