# FASTAFS: toolkit for file system virtualisation of random access compressed FASTA files
----

# Citing FASTAFS or its virtualisation philosophy

[10.1186/s12859-021-04455-3](https://doi.org/10.1186/s12859-021-04455-3)

```
Hoogstrate, Y., Jenster, G.W. & van de Werken, H.J.G.
FASTAFS: file system virtualisation of random access compressed FASTA files.
BMC Bioinformatics 22, 535 (2021).
https://doi.org/10.1186/s12859-021-04455-3
```

----

Direct link to the file format specification:

[https://github.com/yhoogstrate/fastafs/blob/master/doc/FASTAFS-FORMAT-SPECIFICATION.md](https://github.com/yhoogstrate/fastafs/blob/master/doc/FASTAFS-FORMAT-SPECIFICATION.md)

----

Links:
[bio.tools](https://bio.tools/fastafs)

---

![](https://bioinf-galaxian.erasmusmc.nl/public/images/fastafs/fastafs-example.gif)

## Elegant integration of sequence data archives, backwards compatible with FASTA and no API's needed

RNA, DNA and protein sequences are commonly stored in the FASTA format. Although very commonly used and easy to read, FASTA files come with additional metadata files and consume unnecessary disk space. These additional metadata files need to be are necessary to achieve random access and have certain interoperability features, and require additional maintaince. Classical FASTA (de-)compressors only offer back and forwards compression of the files, often requiring to decompress to a new copy of the FASTA file making it inpractical solutions in particular for random access use cases. Although they typically produce very compact archives with quick algorithms, they are not widely adopted in our bioinformatics software.

Here we propose a solution; a virtual layer between (random access) FASTA archives and read-only access to FASTA files and their guarenteed in-sync FAI, DICT and 2BIT files, through the File System in Userspace (FUSE) file system layer. When the archive is mounted, fastafs virtualizes a folder containing the FASTA and necessary metadata files, only accessing the chunks of the archive needed to deliver to the file request. This elegant software solution offers several advantages:
 - virtual files and their system calls are identical to flat files and preserve backwards compatibility with tools only compatible with FASTA, also for random access use-cases,
 - there is no need to use additional disk space for temporary decompression or to put entire FASTA files into memory,
 - for random access requests, computational resources are only spent on decompressing the region of interest,
 - it does not need multiple implementations of software libraries for each distinct tool and for each programming language,
 - it does not require to maintain multiple files that all together make up one data entity as it is guaranteed to provide dict- and fai-files that are in sync with their FASTA of origin.

In addition, the corresponding toolkit offers an interface that allows ENA sequence identification, file integrity verification and  management of the mounted files and process ids.

FASTAFS is deliberately made backwards compatible with both TwoBit and Fasta. The package even allows to mount TwoBit files instead of FASTAFS files, to FASTA files. For those who believe FASTAFS is this famous 15th standard (<https://xkcd.com/927/>)?
Partially, it is not designed to replace FASTA nor TwoBit as the mountpoints provide an exact identical way of file access as regular flat file acces, and is thus backwards compatible. Instead, it offers the same old standard with an elegant toolkit that allows easier integration with workflow management systems.

## Installation and compilation

Currently the package uses cmake for compilation
Required dependencies are:

 -   libboost (only for unit testing, will be come an optional dependency soon)
 -   libopenssl (for generating MD5 hashes)
 -   libfuse (for access to the fuse layer system and file virtualization)
 -   c++ compiler supporting c++-14
 -   glibc
 -   libssl (for checking sequences with ENA)
 -   zlib (crc32 checksum)
 -   cmake or meson + ninja-build
 -   libcrypto for MD5sums
 



```
# debian + ubuntu like systems:

sudo apt install git build-essential cmake libboost-dev libssl-dev libboost-test-dev libboost-system-dev libboost-filesystem-dev zlib1g-dev libzstd-dev libfuse-dev
git clone https://github.com/yhoogstrate/fastafs.git
cd fastafs
```

```
# RHEL + CentOS + Fedora like systems:

sudo yum install git cmake gcc-c++ boost-devel openssl-devel libzstd-devel zlib-devel fuse-devel
git clone https://github.com/yhoogstrate/fastafs.git
cd fastafs
```



Compile (release, recommended):
```
cmake -DCMAKE_BUILD_TYPE=release -DCMAKE_INSTALL_PREFIX=/usr/local -DCMAKE_VERBOSE_MAKEFILE:BOOL=ON .
make "$@" -j `nproc`
sudo make install
```

If you do not have root permission, use the following instead:
```
cmake -DCMAKE_BUILD_TYPE=release -DCMAKE_INSTALL_PREFIX=~/.local -DCMAKE_VERBOSE_MAKEFILE:BOOL=ON .
make "$@" -j `nproc`
make install
```

If you like to play with the code and like to help development, you can create a debug binary as follows:
```
cmake -DCMAKE_BUILD_TYPE=debug -DCMAKE_INSTALL_PREFIX=/usr/local -DCMAKE_VERBOSE_MAKEFILE:BOOL=ON .
make "$@" -j `nproc`
sudo make install
```


If you have patches, changes or even cool new features you believe are worth contributing, please run astyle with the following command:

```
$ make tidy
```

This styles the code in a more or less compatible way with the rest of the code.
Thanks in advance(!)

## usage
### fastafs cache: adding files to fastafs
We can add files to the fastafs database by running:

```
$ fastafs cache test ./test.fa
```

Or, starting with 2bit:

```
$ fastafs cache test-from-2bit ./test.2bit
```

FASTAFS files will be saved in `~/.local/share/fastafs/<uid>.fastafs` and an entry will be added to the 'database' (`~/.local/share/fastafs/index`).

### fastafs list: overview of fastafs db

The `list` command lists all FASTAFS files located in the 'database' (`~/.local/share/fastafs/index`):

```
$ fastafs list

FASTAFS NAME    FASTAFS        SEQUENCES    BASES   DISK SIZE
test            v0-x32-2bit    7            88      214      
```

### fastafs info: stats of cached fasta file
```
$ fastafs info

# FASTAFS NAME: /home/youri/.local/share/fastafs/test.fastafs
# SEQUENCES:    7
chr1                    16          2bit    75255c6d90778999ad3643a2e69d4344
chr2                    16          2bit    8b5673724a9965c29a1d76fe7031ac8a
chr3.1                  13          2bit    61deba32ec4c3576e3998fa2d4b87288
chr3.2                  14          2bit    99b90560f23c1bda2871a6c93fd6a240
chr3.3                  15          2bit    3625afdfbeb43765b85f612e0acb4739
chr4                    8           2bit    bd8c080ed25ba8a454d9434cb8d14a68
chr5                    6           2bit    980ef3a1cd80afec959dcf852d026246
```

### fastafs mount: mount fastafs archive to unlock fasta file(s)

```
$ fastafs mount hg19 /mnt/fastafs/hg19 
$ ls /mnt/fastafs/hg19 
hg19.2bit  hg19.dict  hg19.fa  hg19.fa.fai  seq

$ ls -alsh /mnt/fastafs/hg19
total 0
-rw-r--r-- 1 youri youri 779M Aug 19 15:26 hg19.2bit
-rw-r--r-- 1 youri youri 7.9K Aug 19 15:26 hg19.dict
-rw-r--r-- 1 youri youri 3.0G Aug 19 15:26 hg19.fa
-rw-r--r-- 1 youri youri 3.5K Aug 19 15:26 hg19.fa.fai
drwxr-xr-x 1 youri youri    0 Aug 19 15:26 seq


$ head -n 5 /mnt/bio/hg19/hg19.dict 
@HD	VN:1.0	SO:unsorted
@SQ	SN:chr1	LN:249250621	M5:1b22b98cdeb4a9304cb5d48026a85128	UR:fastafs:///hg19
@SQ	SN:chr2	LN:243199373	M5:a0d9851da00400dec1098a9255ac712e	UR:fastafs:///hg19
@SQ	SN:chr3	LN:198022430	M5:641e4338fa8d52a5b781bd2a2c08d3c3	UR:fastafs:///hg19
@SQ	SN:chr4	LN:191154276	M5:23dccd106897542ad87d2765d28a19a1	UR:fastafs:///hg19
```

### fastafs mount: use custom padding

```
$ fastafs mount test /mnt/fastafs/test        
$ ls /mnt/fastafs/test 
test.2bit  test.dict  test.fa  test.fa.fai

$ cat /mnt/fastafs/test/test.fa
>chr1
ttttccccaaaagggg
>chr2
ACTGACTGnnnnACTG
>chr3.1
ACTGACTGaaaac
>chr3.2
ACTGACTGaaaacc
>chr3.3
ACTGACTGaaaaccc
>chr4
ACTGnnnn
>chr5
nnACTG

$ umount /mnt/fastafs/test

$ fastafs mount -p 4 test /mnt/fastafs/test        
$ cat  /mnt/fastafs/test/test.fa | head -n 15
>chr1
tttt
cccc
aaaa
gggg
>chr2
ACTG
ACTG
nnnn
ACTG
>chr3.1
ACTG
ACTG
aaaa
c
```

To find the file size of chrM (16571):
```
$ ls -l /mnt/bio/hg19/seq/chrM

-rw-r--r-- 1 youri youri 16571 Feb  1 10:47 /mnt/bio/hg19/seq/chrM
```

### Find all running `fastafs mount` / `mount.fastafs` instances

The output format of `fastafs ps` is: `<pid>\t<source>\t<destination>\n`

```
$ fastafs ps
16383	/home/youri/.local/share/fastafs/test.fastafs	/mnt/tmp
```

### Mounting via fstab (for instance on linux boot)

You can add the following line(s) to /etc/fstab to make fastafs mount during boot:

```
mount.fastafs#/home/youri/.local/share/fastafs/hg19.fastafs /mnt/fastafs/hg19 fuse auto,allow_other 0 0
```

Here `mount.fastafs` refers to the binary that only does mounting, without the rest of the toolkit.
This is followed by a hash-tag and the path of the desired fastafs file. The next value is the mount point followed by the argument indicating it runs fuse.
The `auto,allow_other` refers to the `-o` argument.
Here `auto` ensures it is mounted automatically after boot.
Given that a system mounts as root user, the default permission of the mount point will be only for root. 
By setting `allow_other`, file system users get permissions to the mountpoint.

## Contributing
Feel free to start a discussion or to contribute to the GPL licensed code.
If you are willing to make even the smallest contribution to this project in any way, really, feel free to open an issue or to send an e-mail.

[![Codacy Badge](https://api.codacy.com/project/badge/Grade/c90c7d61651d4e18aa82a4b02f3599fa)](https://www.codacy.com/app/yhoogstrate/fastafs?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=yhoogstrate/fastafs&amp;utm_campaign=Badge_Grade)
This should be identical to the reference implementation of zstd-seekable but with minor tweaks to get it working with CPP
# FASTAFS Format specification

The file format consists of four parts:
 - Header, for initiation and format recognition
 - (Compressed) Sequence data
 - Index, metadata required to read sequence data
 - (Optional) metadata

Although it is uncommon, the FASTAFS format denotes the index (a typical header) after the data.
The reason for this location is that most metadata can only be known after converting the entire FASTA file to FASTAFS.
For instance, the number of nucleotides encoded within 2bit data, the number of N-characters or the checksum.
If this metadata would be written in the header located before the sequence data, either:
 - A lot of memory has to be allocated to the heap
 - The FASTA files have to be read twice (once for parsing, once for 2bit conversion)
 - Separate tmp files need to be written to disk and merged afterwards (may cost loads of IO traffic)

## Layout ## 

| Section | Sub | Size | Description |
| ------ | ------ | ------ | ------ |
| GENERIC-HEADER |        |        |        |
|        | [MAGIC](#magic) | 4 bytes | `x0F x0A x46 x53` |
|        | [FILE FORMAT VERSION](#file-format-version) | [4-byte integer](#four-byte-integer) | `x00 x00 x00 x00` |
|        | [FASTAFS-FLAGS](#fastafs-flags) | 2 bytes | Certain binary flags |
|        | [FILE-POSITION-OF-INDEX](#file-position-of-the-index) | [4-byte integer](#four-byte-integer) | Location in the file (offset in bytes from beginning) where the INDEX is located | 
| DATA | --- | --- | --- |
|   -> per sequence | 
|        | N-COMPRESSED-NUCLEOTIDES | uint32_t as [4-byte integer](#four-byte-integer) | The number of 2bit compressed nucleotides. Technical limit is 256^4 |
|        | [TWOBIT-DATA](#twobit-data) | sequence of 2bit-bytes | length can be deduced from header |
|        | UNKNOWN-NUCLEOTIDES | uint32_t as [4-byte integer](#four-byte-integer) | Number of N-entries |
|        | N-STARTS | N x uint32_t as [4-byte integer](#four-byte-integer) | start positions (0-based) |
|        | N-ENDS | N x uint32_t as [4-byte integer](#four-byte-integer) | end positions (0-based) |
|        | [MD5-CHECKSUM](#md5-checksum) | 16 x byte | MD5 compatible with CRAM, BAM, DICT & ENA |
|        | ** RESERVED-REGIONS | uint32_t as [4-byte integer](#four-byte-integer) | Number of R-entries (reserved regions ~ incomplete file) - not yet implemented and must be enabled by a flag |
|        | ** R-STARTS | N x uint32_t as [4-byte integer](#four-byte-integer) | start positions (0-based) |
|        | ** R-ENDS | N x uint32_t as [4-byte integer](#four-byte-integer) | end positions (1-based) |
|        | MASKED-NUCLEOTIDES | uint32_t as [4-byte integer](#four-byte-integer) | Number of M-entries (for lower case regions) |
|        | M-STARTS | M x uint32_t as [4-byte integer](#four-byte-integer) | start positions (0-based) - default is CAPITAL, within M-blocks is LOWER-case |
|        | M-ENDS | M x uint32_t as [4-byte integer](#four-byte-integer) | end positions (0-based) |
| INDEX  | --- | --- |  |
|        | NUMBER-SEQUENCES | uint32_t as [4-byte integer](#four-byte-integer) | Number of sequences included |
|   -> per sequence | 
|        | [SEQUENCE-FLAGS](#sequence-flags) | 2 bytes | storing metadata and type of data |
|        | NAME-LENGTH | 1 byte as unsigned char | length in bytes; name cannot exceed 255 bytes |
|        | NAME-FASTA | NAME-LENGTH x char | FASTA header; may not include new-lines or '>' |
|        | START-POSITION-IN-BODY of N-COMPR-NUC | uint32_t as [4-byte integer](#four-byte-integer) | Location in the file (offset in bytes from beginning) where the DATA block for this sequence starts |
| METADATA | only first char is required, the rest is always optional |
|          | N-METADATA-TAGS | 1 x char |
| METADATA-ENTRY [per entry] |  ~ limits to 'only' 256 distinct types of metadata
|          | METADATA-TYPE-FLAG | 2 bytes | 
|          | ENTRY | type specific, examples below: |
|          | => ORIGINAL PADDING | uint32_t as [4-byte integer](#four-byte-integer) | The number of nucleotides per line in the original FASTA file |
| CRC32  | Checksum on entire file | 4 bytes | To ensure whole file integrity |

### GENERIC-HEADER ###

#### Magic ####

The file magic (first four bytes used to recognise binary file types) of FASTAFS are `x0F x0A x46 x53`.
They stand for `FA` (as HEX) and `FS` (in ASCII).
The bit representation of these bytes are:

```
    +--------+--------+--------+--------+
    |00001111|00001010|01000110|01010011|
    +--------+--------+--------+--------+
```

#### FILE FORMAT VERSION ####

The version of the file format specification implemented as four byte integer. Currently only one version exists: 

`x00 x00 x00 x00`

The bit representation of these bytes are:

```
    +--------+--------+--------+--------+
    |00000000|00000000|00000000|00000000|
    +--------+--------+--------+--------+
```

#### FASTAFS-FLAGS ####

``` 
bit 0   file-complete
bit 1   reserved [most likely for 1 = 64-bit fastafs]
bit 2   reserved
bit 3   reserved
bit 4   reserved
bit 5   reserved
bit 6   reserved
bit 7   reserved
---
bit 8   reserved
bit 9   reserved
bit 10  reserved
bit 11  reserved
bit 12  reserved
bit 13  reserved
bit 14  reserved
bit 15  reserved
```

File-complete set to 1 means that the file writing has completed.
If this value is set to 0, the file is either being written or corrupt (because interrupted write process)


#### File position of the index ####

The index is located at the end of the data. This file offset in bytes from the files start position is indicated as [4-byte integer](#four-byte-integer).


### DATA ###

Repeated for every sequence, in order matching SEQUENCE-HEADER

#### SEQUENCE-FLAGS #### 

The sequence flag allows to describe the following metadata for each sequence:

```
bit 0   combined sequence type
bit 1   combined sequence type
```

| bit-0 | bit-1 | Type | Alphabet |
| ---- | ---- | - | - |
| `0` | `0` | DNA | `ACTG` + `N` |
| `1` | `0` | RNA | `ACUG` + `N` |
| `0` | `1` | IUPEC Nucleotide | `ACGTURYKMSWBDHVN` + `-` |
| `1` | `1` | reserved for protein | to be determined |

```
bit 2   reserved    [reserved, library type 2 -> protein]
bit 3   is-complete [1: checksum is present, 0: some regions are reserved but not yet 'downloaded']
bit 4   is-circular 
bit 5   reserved    [reserved for vertical duplicate? located in other file]
bit 6   reserved
bit 7   reserved
---
bit 8   reserved
bit 9   reserved
bit 10  reserved
bit 11  reserved
bit 12  reserved
bit 13  reserved
bit 14  reserved
bit 15  reserved
```

#### TwoBit Data ####

The TwoBit data is encoded in a long array of bytes in which each byte encodes four nucleotides.

The following bits to nucleotide encoding is used:

| bits | Nucleotide |
| ---- | - |
| `00` | T |
| `01` | C |
| `10` | A |
| `11` | G |

Encoded into a byte in the following order:

```
    +----------+------+
    | 00000011 | TTTG |
    +----------+------+
    | 00001100 | TTGT |
    +----------+------+
    | 00110000 | TGTT |
    +----------+------+
    | 11000000 | GTTT |
    +----------+------+
```

#### MD5 checksum ####

Per sequence, an MD5 checksum is stored as it's binary encoded digest.

The MD5 checksum is calculated as described in the CRAM specification (section 11): https://samtools.github.io/hts-specs/CRAMv3.pdf

Only nucleotides are included (no newlines, no whitespaces, no headers), all in uppercase.


### INDEX ###

The index is put to the end because lots of information is unknown during conversion and requires large amounts of RAM allcoated

#### Four Byte Integer ####

A four byte integer is a binary encoded integer value (using 4 bytes).

 - A 0 is encoded as follows (`x00 x00 x00 x00`):

```
    +--------+--------+--------+--------+
    |00000000|00000000|00000000|00000000|
    +--------+--------+--------+--------+
```

 - A 1 is encoded as follows (`x00 x00 x00 x01`):

```
    +--------+--------+--------+--------+
    |00000000|00000000|00000000|00000001|
    +--------+--------+--------+--------+
```

 - A 2 is encoded as follows (`x00 x00 x00 x02`):

```
    +--------+--------+--------+--------+
    |00000000|00000000|00000000|00000010|
    +--------+--------+--------+--------+
```

 - A 3 is encoded as follows (`x00 x00 x00 x03`):

```
    +--------+--------+--------+--------+
    |00000000|00000000|00000000|00000011|
    +--------+--------+--------+--------+
```

 - A 255 is encoded as follows (`x00 x00 x00 xFF`):

```
    +--------+--------+--------+--------+
    |00000000|00000000|00000000|11111111|
    +--------+--------+--------+--------+
```

 - A 256 is encoded as follows (`x00 x00 x01 x00`):

```
    +--------+--------+--------+--------+
    |00000000|00000000|00000001|00000000|
    +--------+--------+--------+--------+
```

 - A 257 is encoded as follows (`x00 x00 x01 x01`):

```
    +--------+--------+--------+--------+
    |00000000|00000000|00000001|00000001|
    +--------+--------+--------+--------+
```

Denote that this implementation is DIFFERENT than for the [UCSC implementation](#four-byte-integer-ucscu-implementation)!


#### Four Byte Integer UCSC implementation ####

A four byte integer is a binary encoded integer value (using 4 bytes). 


 - A 0 is encoded as follows (`x00 x00 x00 x00`):

```
    +--------+--------+--------+--------+
    |00000000|00000000|00000000|00000000|
    +--------+--------+--------+--------+
```

 - A 1 is encoded as follows (`x01 x00 x00 x00`) - DENOTE INVERSED ORDER OF INTEGERS:

```
    +--------+--------+--------+--------+
    |00000001|00000000|00000000|00000000|
    +--------+--------+--------+--------+
```

 - A 2 is encoded as follows (`x02 x00 x00 x00`) - DENOTE INVERSED ORDER OF INTEGERS:

```
    +--------+--------+--------+--------+
    |00000010|00000000|00000000|00000000|
    +--------+--------+--------+--------+
```

 - A 3 is encoded as follows (`x03 x00 x00 x00`) - DENOTE INVERSED ORDER OF INTEGERS:

```
    +--------+--------+--------+--------+
    |00000011|00000000|00000000|00000000|
    +--------+--------+--------+--------+
```

 - A 255 is encoded as follows (`xFF x00 x00 x00`) - DENOTE INVERSED ORDER OF INTEGERS:

```
    +--------+--------+--------+--------+
    |11111111|00000000|00000000|00000000|
    +--------+--------+--------+--------+
```

 - A 256 is encoded as follows (`x00 x01 x00 x00`) - DENOTE INVERSED ORDER OF INTEGERS:

```
    +--------+--------+--------+--------+
    |00000000|00000001|00000000|00000000|
    +--------+--------+--------+--------+
```

 - A 257 is encoded as follows (`x01 x01 x00 x00`) - DENOTE INVERSED ORDER OF INTEGERS:

```
    +--------+--------+--------+--------+
    |00000001|00000001|00000000|00000000|
    +--------+--------+--------+--------+
```


Denote that this implementation is DIFFERENT than for the FASTAFS implementation!
 

# TODO's #
Add SEQUENCE-TYPE (1 byte):
  * ACTG / ACUG
  * sequence size:
    - 256^4 (uint)  <- current default and only option
    - 256^8 (samtools (cram) ITF8)

Add DEFAULT-PADDING to GENERIC-HEADER

Add tqdm like progress - if possible
