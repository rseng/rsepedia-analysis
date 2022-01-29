[![Build Status](https://travis-ci.com/herrsalmi/FConverter.svg?branch=master)](https://github.com/herrsalmi/FConverter/releases/latest)
[![codecov](https://codecov.io/gh/herrsalmi/FConverter/branch/master/graph/badge.svg)](https://codecov.io/gh/herrsalmi/FConverter)
# FConverter
File format converter for FImpute2

### Running FConverter
```
usage : FConverter <vcf2fimpute|snpID|fimpute2vcf> [options]
	help | -h
		print this help and exit the program
	-v
		print the program version
	vcf2fimpute
		vcf=[file]          gzip compressed vcf file 
		chip=[integer]      chip number
		[-p]                print file header
		[nthr=xx]           use xx threads (default=4)
	snpID
		hd=[file]           gzip compressed vcf file
		[ld=[file]]         gzip compressed vcf file (comma separated if many)
	fimpute2vcf
		gen=[file]          imputed genotypes file 
		snp=[file]          snp info file reported by FImpute 
		vcf=[file]          gzip compressed vcf file (reference) (comma separated if many)
		out=[string]        output file prefix
		chip=[integer]      chip number (0 for a VCF with all chips)

```
