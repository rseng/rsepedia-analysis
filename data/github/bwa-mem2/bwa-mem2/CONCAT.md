[![GitHub Downloads](https://img.shields.io/github/downloads/bwa-mem2/bwa-mem2/total?label=GitHub%20Downloads)](https://github.com/bwa-mem2/bwa-mem2/releases)
[![BioConda Install](https://img.shields.io/conda/dn/bioconda/bwa-mem2?label=BioConda%20Installs)](https://anaconda.org/bioconda/bwa-mem2)

## Important Information

***We are happy to announce that the index size on disk is down by 8 times and in memory by 4 times due to moving to only one type of FM-index (2bit.64 instead of 2bit.64 and 8bit.32) and 8x compression of suffix array. For example, for human genome, index size on disk is down to ~10GB from ~80GB and memory footprint is down to ~10GB from ~40GB.***
***There is a substantial reduction in index IO time due to the reduction and hardly any performance impact on read mapping.***
***Due to this change in index structure (in commit #4b59796, 10th October 2020), you will need to rebuild the index.***

***Added MC flag in the output sam file in commit a591e22. Output should match original bwa-mem version 0.7.17.***

***As of commit e0ac59e, we have a git submodule safestringlib. To get it, use --recursive while cloning or use "git submodule init" and "git submodule update" in an already cloned repository (See below for more details).***


## Getting Started
```sh
# Use precompiled binaries (recommended)
curl -L https://github.com/bwa-mem2/bwa-mem2/releases/download/v2.0pre2/bwa-mem2-2.0pre2_x64-linux.tar.bz2 \
  | tar jxf -
bwa-mem2-2.0pre2_x64-linux/bwa-mem2 index ref.fa
bwa-mem2-2.0pre2_x64-linux/bwa-mem2 mem ref.fa read1.fq read2.fq > out.sam

# Compile from source (not recommended for general users)
# Get the source
git clone --recursive https://github.com/bwa-mem2/bwa-mem2
cd bwa-mem2
# Or
git clone https://github.com/bwa-mem2/bwa-mem2
cd bwa-mem2
git submodule init
git submodule update
# Compile and run
make
./bwa-mem2
```

## Introduction

Bwa-mem2 is the next version of the bwa-mem algorithm in [bwa][bwa]. It
produces alignment identical to bwa and is ~1.3-3.1x faster depending on the use-case, dataset and the running machine.

The original bwa was developed by Heng Li (@lh3). Performance enhancement in
bwa-mem2 was primarily done by Vasimuddin Md (@yuk12) and Sanchit Misra (@sanchit-misra)
from Parallel Computing Lab, Intel.
Bwa-mem2 is distributed under the MIT license.

## Installation

For general users, it is recommended to use the precompiled binaries from the
[release page][rel]. These binaries were compiled with the Intel compiler and
runs faster than gcc-compiled binaries. The precompiled binaries also
indirectly support CPU dispatch. The `bwa-mem2` binary can automatically choose
the most efficient implementation based on the SIMD instruction set available
on the running machine. Precompiled binaries were generated on a CentOS7
machine using the following command line:
```sh
make CXX=icpc multi
```

[bwa]: https://github.com/lh3/bwa
[rel]: https://github.com/bwa-mem2/bwa-mem2/releases

## Usage

The usage is exactly same as the original BWA MEM tool. Here is a brief synopsys. Run ./bwa-mem2 for available commands.

```sh
# Indexing the reference sequence (Requires 28N GB memory where N is the size of the reference sequence).
./bwa-mem2 index [-p prefix] <in.fasta>
Where 
<in.fasta> is the path to reference sequence fasta file and 
<prefix> is the prefix of the names of the files that store the resultant index. Default is in.fasta.

# Mapping 
# Run "./bwa-mem2 mem" to get all options
./bwa-mem2 mem -t <num_threads> <prefix> <reads.fq/fa> > out.sam
Where <prefix> is the prefix specified when creating the index or the path to the reference fasta file in case no prefix was provided.
```

## Performance

Datasets:  
Reference Genome: human_g1k_v37.fasta

 Alias	    |  Dataset source				|  No. of reads	| Read length 
 --------- | --------- | --------- | --------- 
 D1	|  Broad Institute				|  2 x 2.5M	bp	|	151bp
 D2	|  SRA: SRR7733443				|  2 x 2.5M	bp	|	151bp  
 D3	|  SRA: SRR9932168				|  2 x 2.5M	bp	|	151bp  
 D4	|  SRA: SRX6999918				|  2 x 2.5M	bp	|	151bp  



Machine details:  
Processor: Intel(R) Xeon(R) 8280 CPU @ 2.70GHz  
OS: CentOS Linux release 7.6.1810  
Memory: 100GB  


We followed the steps below to collect the performance results:  
A. Data download steps:
1. Download SRA toolkit from https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=software#header-global    
2. tar xfzv sratoolkit.2.10.5-centos_linux64.tar.gz  
3. Download D2: sratoolkit.2.10.5-centos_linux64/bin/fastq-dump --split-files SRR7733443   
4. Download D3: sratoolkit.2.10.5-centos_linux64/bin/fastq-dump --split-files SRR9932168   
5. Download D4: sratoolkit.2.10.5-centos_linux64/bin/fastq-dump --split-files SRX6999918   



B. Alignment steps:   
1. git clone https://github.com/bwa-mem2/bwa-mem2.git   
2. cd bwa-mem2   
3. ```make CXX=icpc``` (using intel C/C++ compiler)   
or   ```make``` (using gcc compiler)   
4. ./bwa-mem2 index <ref.fa>   
5. ./bwa-mem2 mem [-t <#threads>] <ref.fa> <in_1.fastq> [<in_2.fastq>]  >  <output.sam>   

For example,  in our double socket (56 threads each) and double numa compute node, we used the following command line to align D2 to human_g1k_v37.fasta reference genome.  
```
numactl -m 0 -C 0-27,56-83 ./bwa-mem2 index human_g1k_v37.fasta  
numactl -m 0 -C 0-27,56-83 ./bwa-mem2 mem -t 56 human_g1k_v37.fasta SRR7733443_1.fastq SRR7733443_2.fastq > d3_align.sam
```

<p align="center">
<img src="https://github.com/bwa-mem2/bwa-mem2/blob/master/images/bwa-mem2-1.png" height="400"/a></br>
<img src="https://github.com/bwa-mem2/bwa-mem2/blob/master/images/bwa-mem2-2.png" height="400"/a></br>
<img src="https://github.com/bwa-mem2/bwa-mem2/blob/master/images/bwa-mem2-3.png" height="400"/a></br>
<img src="https://github.com/bwa-mem2/bwa-mem2/blob/master/images/bwa-mem2-4.png" height="400"/a></br>
</p> 

## bwa-mem2 seeding speedup with Enumerated Radix Trees (Code in ert branch)

The ert branch of bwa-mem2 repository contains codebase of enuerated radix tree based acceleration of bwa-mem2. The ert code is built on the top of bwa-mem2 (thanks to the hard work by @arun-sub). 
The following are the highlights of the ert based bwa-mem2 tool: 
1. Exact same output as bwa-mem(2) 
2. The tool has two additional flags to enable the use of ert solution (for index creation and mapping), else it runs in vanilla bwa-mem2 mode 
3. It uses 1 additional flag to create ert index (different from bwa-mem2 index) and 1 additional flag for using that ert index (please see the readme of ert branch) 
4. The ert solution is 10% - 30% faster (tested on above machine configuration) in comparison to vanilla bwa-mem2 -- users are adviced to use option `-K 1000000` to see the speedups 
5. The memory foot print of the ert index is ~60GB 
6. The code is present in ert branch: https://github.com/bwa-mem2/bwa-mem2/tree/ert


## Citation

Vasimuddin Md, Sanchit Misra, Heng Li, Srinivas Aluru.
<b> Efficient Architecture-Aware Acceleration of BWA-MEM for Multicore Systems. </b>
<i> IEEE Parallel and Distributed Processing Symposium (IPDPS), 2019. </i>

# bwa-mem2
----------------------------------------------------------------------------------
IMPORTANT:  
1.'bwa-mem2 index' requires ~48GB of space for full human genome plus index.  
2. We recommend to run the code in AVX512 mode to get maximum performance.  
3. The command-line of bwa-mem2 is exactly same as bwa-mem.  
----------------------------------------------------------------------------------


Compile:  
$ ```make```  
'make' compiles using SSE4.1 vector mode by default.

Compile options:  
```make CXX=<compiler>```, e.g compiler can 'g++' or 'icpc'  
```make arch=<mode>```  

mode:  
a. native - generates 'bwa-mem2' binary for the architechture on which it is compiled. It detects the underlying harware flags for AVX512/AVX2/SSE2 vector modes and compiles accordingly. 
For example, Intel Xeon Haswell supports AVX2 vector mode; Intel Xeon Skylake supports AVX512 vector mode. AVX vector support by the processor can be checked in /proc/cpuinfo file.  
b. sse - generates 'bwa-mem2' SSE4.1 vector code binary  
c. avx2 - generates 'bwa-mem2' AVX2 vector code binary  
d. avx512bw - generates 'bwa-mem2' AVX512BW vector code binary  


 ```make multi```  
Generates binaries for all the three vector modes: bwa-mem2.sse41, bwa-mem2.avx2, bwa-mem2.avx512bw and also generates a binary called bwa-mem2 that, if run, detects the underlying architecture and runs the corresponding binary out of bwa-mem2.sse41, bwa-mem2.avx2, bwa-mem2.avx512bw.  


Run:   
   - Create the index: 'bwa-mem2 index <reference file>'  
   - Mapping: 'bwa-mem2 mem' can be run with appropriate paramters  
   e.g. bwa-mem2 mem -t 1 <index> <read1.fq> [<read2.fq>] > <out.sam>  

Notes:  
- In SSE4.1 vector mode, `bwa-mem2` uses bwt index with 2-bit representation  
- In AVX512/AVX2 vector mode, `bwa-mem2` uses bwt index with 8-bit representation  


Example run:  
```sh bwa-mem2 index datasets/hg38.fa```   
```bwa-mem2 mem -t 1 datasets/hg38.fa datasets/SRR099967.filt.fastq > datasets/aln-se.sa```



NUMA:  
    In multi-socket system, it is beneficial to bind the process memory to a socket  
    Linux command `numactl` can be used to bind process memory to a given numa domain  
    For exmaple, the following command binds the memory to socket-0  
``` $ numactl -m 0 bwa-mem2 mem -t 10 -o datasets/aln-se.sa datasets/hg38.fa datasets/SRR099967.filt.fastq```  

Similarly, the process threads can be bound to cores.  
```$ numactl -m 0 -C 0-9 bwa-mem2 mem -t 10 -o datasets/aln-se.sa datasets/hg38.fa datasets/SRR099967.filt.fastq```  

```$ lscpu``` - Linux command provides the architectural details such as list of NUMA domains and their CPU lists  


Directory structure (bwa-mem2):  
- src/ - contains source code files  
- test/ - contains test files to run individual benchmarks for smem, sal, bsw (disabled for now)  
- ext/ - contains external Intel safestring library for secure string operations  
- images/ - contains bwa-mem2 perfomance charts, demonstrating its speedup against original bwa-mem  


Notes:  
- As mentioned previously, the code supports four exceution modes:  
  1. Vector SSE2: It runs bwa-mem2 with the kernels running SSE2 vector instructions (128 bits vector-width)  
  2. Vector AVX2: It runs bwa-mem2 with the kernels running AVX2 vector instructions (256 bits vector-width)  
  3. Vector AVX512: It runs bwa-mem2 with the kernels running AVX512 vector instructions (512 bits vector-width)  

- BWA-MEM2 displays run-time performance profiling of the code on standard output. (Hopefully, we managed it to be self-explanatory)  

----------------------------------------------------------------------------------
Release 2.2.1 (17 March 2021)
---------------------------------
Hotfix for v2.2: Fixed the bug mentioned in #135.


Release 2.2 (8 March 2021)
---------------------------------
Changes since the last release (2.1):

* Passed the validation test on ~88 billions reads (Credits: Keiran Raine, CASM division, Sanger Institute)
* Fixed bugs reported in #109 causing mismatch between bwa-mem and bwa-mem2
* Fixed the issue (# 112) causing crash due to corrupted thread id 
* Using all the SSE flags to create optimized SSE41 and SSE42 binaries


Release 2.1 (16 October 2020)
---------------------------------
Release 2.1 of BWA-MEM2.

Changes since the last release (2.0):
* *Smaller index*: the index size on disk is down by 8 times and in memory by 4 times due to moving to only one type of FM-index (2bit.64 instead of 2bit.64 and 8bit.32) and 8x compression of suffix array. For example, for human genome, index size on disk is down to ~10GB from ~80GB and memory footprint is down to ~10GB from ~40GB. There is a substantial decrease in index IO time due to the reduction and hardly any performance impact on read mapping.

* Added support for 2 more execution modes: sse4.2 and avx.

* Fixed multiple bugs including those reported in Issues #71, #80 and #85.

* Merged multiple pull requests.


Release 2.0 (9 July 2020)
---------------------------------
This is the first production release of BWA-MEM2.

Changes since the last release:
* Made the source code more secure with more than 300 changes all across it.

* Added support for memory re-allocations in case the pre-allocated fixed memory is insufficient.

* Added support for MC flag in the sam file and support for -5, -q flags in the command line.

* The output is now identical to the output of bwa-mem-0.7.17.

* Merged index building code with FMI_Search class.

* Added support for different ways to input read files, now, it is same as bwa-mem.

* Fixed a bug in AVX512 sam processing part, which was leading to incorrect output.


Release 2.0pre2 (4 February 2020)
---------------------------------

Miscellaneous changes:

 * Changed the license from GPL to MIT.

 * IMPORTANT: the index structure has changed since commit 6743183. Please
   rebuild the index if you are using a later commit or the new release.

 * Added charts in README.md comparing the performance of bwa-mem2 with bwa-mem.

Major code changes:

 * Fixed working for variable length reads.

 * Fixed a bug involving reads of length greater than 250bp.

 * Added support for allocation of more memory in small chunks if large
   pre-allocated fixed memory is insufficient. This is needed very rarely
   (thus, having no impact on performance) but prevents asserts from failing
   (code from crashing) in that scenario. 

 * Fixed a memory leak due to not releasing the memory allocated for seeds
   after smem.

 * Fixed a segfault due to non-alignment of small allocated memory in the
   optimized banded Smith-Waterman.

 * Enabled working with genomes larger than 7-8 billion nucleotides (e.g. Wheat
   genome).

 * Fixed a segfault occuring (with gcc compiler) while reading the index.
