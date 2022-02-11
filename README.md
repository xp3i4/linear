Linear 
====
![example workflow](https://github.com/xp3i4/linear/actions/workflows/cmake.yml/badge.svg)
[![License](https://img.shields.io/badge/License-BSD%203--Clause-blue.svg)](https://opensource.org/licenses/BSD-3-Clause)
![platforms](https://img.shields.io/badge/platform-linux-informational.svg)

## ALIgNment-freE method for long-read vARiants resolution 
Structural variants (SVs) detection pipelines for 3Gen sequencing are commonly based on the alignment.
However, rigid static gap models in alignment are less flexible and efficient for diverse types of SVs.
Thus we develop the alignment-free method for SVs detection.
The alignment-free results can be called directly by alignment based tools such as the SVs caller PBSV.

## Build and usage 
### Prerequisites
Please make sure the following systems have been installed before building from the source.
#### Linux only:
||requirement|version|comment|
|--|--|--|--|
|**Compiler**|[GCC](https://gcc.gnu.org/) | ≥ 4.9.0|no other compiler is currently supported|
|**Build system**| [CMAKE](https://cmake.org/) |≥ 3.0.0|/|
|**External libs**|[zlib](https://github.com/madler/zlib)|≥ 1.2|required for `*.gz` and `*.bam` file support |

To install zlib on
Debian-based: Ubuntu, etc
```bash
sudo apt-get install zlib1g zlib1g-dev
```
RedHat-based: Fedora, etc
```bash
sudo dnf install zlib-devel
```
#### Windows and osx are currently not supported.
### Build
To build from the source, create a new directory. In the new directory:
```bash
$CMake [path to source] 
$make linear -j 8 
```
Note:The option `-j 8` is to set up 8 threads to speedup the compilation.

### Usage
Supported formats  for input: `*.fa(stq)(.gz)`.
```bash
$linear read.fastq genome.fa
``` 
Please add argument <b>'-x'</b> between the reads and references to input more than 2 files. 
```bash
$linear *.fastq x *.fa
``` 
Use -h for more details of options
```bash
$linear -h
```


## Adaption to existing pipelines 🐦
Alignment -free results by Linear currently have been tested with the following pipelines.
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
Alignment-free model for variants is flexible to extend.
We are updating models in Linear continuously if there are better ones.
This can lead to different results between versions.

## Format of alignment-free results 
### SAM/BAM
We extended SAM/BAM for alignment-free results.\
It's compatible with the standard SAM/BAM for alignment. Each field is defined in the following table:
- The 6th column of cigar is redefined in the extended SAM/BAM.
Alignment-free cigar denotes the virtual alignment between 2 points, which is always in the pair of 'MG', where 'M' is 'X' or '=' and 'G' is 'I' and 'D'.
- The 10th column of SEQ is subsequence from read or reference.
- The 12th column of the tag 'SA:Z' is redefined.
Other tags are identical to the standard tag, which can be found at [SAM/BAM format](https://samtools.github.io/hts-specs/SAMv1.pdf) and [Optional tags](https://samtools.github.io/hts-specs/SAMtags.pdf).

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



### Alignment-free mapping file (.APF) 
.APF is a nonstandard format based on the [.PAF](https://github.com/lh3/miniasm/blob/master/PAF.md).
The format is still being improved🛠 to provide more readable alignment-free results.
The file can be enabled/disabled with the option '<b>-ot</b>' in Linear.
The .APF format contains the header(H) and records(R) defined in the following table.
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
|R2|QSTR|Query template mapped base|int|
|R3|RSTR|Reference sequence mapped base|int|
|R4|DY|Distance of R3  to last R3|int|
|R5|DX|Distance of R4  to last R4|int|
|R6|RSTRD|record strand|{'+','-'}|
