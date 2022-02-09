
Linear <img width="60px" src="images/linear_logo-1.svg"/>
====
![example workflow](https://github.com/catx1024/linear/actions/workflows/cmake.yml/badge.svg)
[![License](https://img.shields.io/badge/License-BSD%203--Clause-blue.svg)](https://opensource.org/licenses/BSD-3-Clause)
![platforms](https://img.shields.io/badge/platform-linux-informational.svg)

## ALIgNment-freE method for long-read vARiants resolution 
Structural variants (SVs) pipelines commonly rely on the alignment.
However, rigid static gap model in alignment is inflexible to resolve diverse SVs and inefficient to compute.
Thus we develop the alignment-free method for more efficient dection of SVs in 3rd sequencing reads.
The results can be called directly by alignment based tools such as the SVs caller PBSV and the visualization tool IGV.

## Build and usage 💩
Linear is easy to build and use.
Make sure the following tools or libraries are available before building the source.

- <img src="images/linux_logo.png" width="26"/> GNU/Linux kernel> 4.9.0  GCC > 4.9
- <img src="images/cmake_logo.png" width="26"/> CMAKE > 3.0.0
- <img src="images/Zlib_3D_green.svg" width="36"/> zlib

To install zlib on
Debian-based: Ubuntu, etc
```bash
sudo apt-get install zlib1g zlib1g-dev
```
RedHat-based: Fedora, etc
```bash
sudo dnf install zlib-devel
```
To build Linear from the source create a new directory. In the new directory type:
```bash
$CMake [path to source] 
$make linear -j 8 
```
Note:The <b>'-j 8'</b> option is to set up 8 threads to speedup the compilation.
### Usage
Supported file formats  for input: .fa(.gz) and .fastq(.gz).
```bash
$linear read.fa genome.fa
``` 
Please add argument <b>'-x'</b> if computing more than one reads and genomes.
```bash
$linear *fa x *fa
``` 
Use -h for more details of options
```bash
$linear -h
```

## Format of files 
Output of Linear extends standard SAM/BAM format to alignment-free results.
Definition of each fields are :
### SAM/BAM

|col |filed|Description|Status|
|--|--|--|--|
|   1  | QNAME | Query template NAME                       | Yes       |           
|   2  | FLAG  | bitwise FLAG                              | Yes       | 
|   3  | RNAME | Reference sequence NAME                   | Yes       | 
|   4  | POS   | 1-based leftmost mapping POSition         | Yes       | 
|   5  | MAPQ  | MAPping Quality                           | Yes       | 
|   6  | CIGAR | CIGAR string                              | Redefined   | 
|   7  | RNEXT | Reference name of the mate/next read      | Yes       |
|   8  | PNEXT | Position of the mate/next read            | Yes       |
|   9  | TLEN  | observed Template LENgth                  | Yes       | 
|   10 | SEQ   | segment SEQuence                          | Redefined   |
|   11 | QUAL  | ASCII of Phred-scaled base QUALity+33     | Yes       |
|   12 | TAG   | Optional tags                             | Redefined   |

- The 6th column of cigar is redefined in the extended SAM/BAM.
Alignment-free cigar denotes the virtual alignment between 2 points, which is always in the pair of 'MG', where 'M' is 'X' or '=' and 'G' is 'I' and 'D'.


- The 10th column of SEQ is inferred according to the both read and reference.
The following table lists bases in SEQ from read or reference for each cigar operation, where
1/0 are Yes/No and 0.5 is conditional Yes

|From|=|X|M|I|D|S|H|
|--|--|--|--|--|--|--|--|
|read|0|0.5|0|1|0|1|0|
|ref|1|0|1|0|0|0|0|

- The 12th column is redefined because the tag 'SA:Z' uses cigar.
Other tags are identical to the standard.
Definition of each tag can be found at [SAM/BAM format](https://samtools.github.io/hts-specs/SAMv1.pdf) and [Optional tags](https://samtools.github.io/hts-specs/SAMtags.pdf)

### Approximate mapping file (.apf) 
This is a nonstandard format based on the [.PAF](https://github.com/lh3/miniasm/blob/master/PAF.md)
The format is being developed.🛠
It's to provide readable alignment-free result.
The file can be enabled/disabled with the option '<b>-ot</b>'.
The .apf file contains the header and records.
The following 2 tables lists the fields of the header(H) and record(R).
|col(H/R) |filed|Description|Type|
|--|--|--|--|
|H1|@|sign of header|{'@'}|
|H2| QNAME|Query template NAME|string|
|H3| QLEN|Query template LENGTH|int|
|H4| QSTR|Query template mapped START| int |
|H5| QEND|Query template mapped END| int |
|H5| QSTRD|Query template mapped main STRAND|{'+','-'}|
|H6| RNAME | Reference sequence NAME|String| 
|H6| RLEN | Reference sequence LENGTH|int| 
|H7| RSTR | Reference sequence mapped START|int| 
|H8| REND | Reference sequence mapped END|int| 
|R1|\||sign to start record|string|
|R2|QSTR|Query region start|int|
|R3|RSTR|Reference region start|int|
|R4|DY|Distance of current 3th col to last|int|
|R5|DX|Distance of current 4th col to last|int|
|R6|RSTRD|region strand|{'+','-'}|

Following is an example of the apf the read mapped to the reference.
```
@ S1_14253 10049 45 10049 + chr1 44948 9964 19212 
| 45 9964 0 0 1 - 
| 215 10119 170 155 2 -
| 352 10240 137 121 3 -
| 446 10327 94 87 4 - 
...
| 6656 16080 144 128 48 - 
| 6800 16208 144 128 49 - 
| 6944 16336 144 128 50 - 
| 1931 16507 -5013 171 51 + 
| 2126 16690 195 183 52 + 
| 2309 16857 183 167 53 + 
| 2462 16999 153 142 54 + 
| 2592 17120 130 121 55 + 
| 2674 17189 82 69 56 + 
...
| 9680 18839 0 7 68 - 
| 9870 19020 190 181 69 -
```

## Adaption to existing pipelines 🐦
Currently Linear in SAM has been tested with the following tools or pipelines.
### Adaption to Samtools
The compatibility of alignment results to Samtools 1.10 has been test.
Alignment-free results can work with 'samtools view', 'samtools index' and 'samtools sort' to convert and index SAM/BAM file.

### Adaption to SVs callers
The compatibility of the alignment-free results to the SVs caller [PBSV](https://github.com/PacificBiosciences/pbsv) 2.6.2 has been tested with PacBio raw reads and CCS reads.
SAM/BAM from Linear can work with 'PBSV discover' and 'PBSV call' provided the sample and group name is set appropriately with the -s option in pbsv discover.

### Adaption to seqeunce graphical tools (IGV)
The compatibility of the alignment-free results to the IGV 2.8.3 has been tested.
Alignment-free .bam  can be visualised directly in IGV.

## Updating variant models 🐢
Alignment-free model for variants is flexible to replace.
We are updating models in Linear continuously such as new models for nested variants.
This probably leads to different results between versions.
