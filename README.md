
Linear <img width="80px" src="images/linear_logo-1.svg"/>
====
![example workflow](https://github.com/catx1024/linear/actions/workflows/cmake.yml/badge.svg)
[![License](https://img.shields.io/badge/License-BSD%203--Clause-blue.svg)](https://opensource.org/licenses/BSD-3-Clause)
![platforms](https://img.shields.io/badge/platform-linux-informational.svg)

## ALIgNment-freE method for long-read vARiants resolution 
Structural variants (SVs) pipelines commonly rely on the alignment.
However, rigid static gap model in alignment is inflexible to resolve diverse SVs and inefficient to compute.
Thus we develop the alignment-free method for mor efficient dection of SVs in 3rd sequencing reads.

### Efficient
With diverse algorithms and optimizations, alignment-free methods is commonly orders of magnitude faster with comparable accuracy and higher sensitivity.  

### Flexible
Linear applies more flexible SVs models than these in aligners.
Models includes but not limited to basic types, such as ins, dels, invs, as well as nested ones. 
More enhanced models are flexible to develop and extend.

### Compatible 
The alignment-free results is compatible to the alignment based pipelines.
The results can be called directly by existing tools such as the SVs caller PBSV and the visualization tool IGV.

## Build and usage
Linear is easy to build and use.
Make sure the following tools or libraries are available before building the source.
| Prerequisites  |   Versions          |
| ------------------- | ------------------------- |
|<img src="images/linux_logo.png" width="16"/> Linux| >4.9.0|
|<img src="images/gcc_logo.png" width="16"/> GCC|>4.9|
|<img src="images/cmake_logo.png" width="16"/> CMAKE|>3.0.0|
|<img src="images/Zlib_3D_green.svg" width="16"/> zlib|/|


To build from the source create a new directory. In the directory type:
```bash
$CMake [path to source] 
$make linear 
```
### Usage
Supported file formats  for input: .fa(.gz) and .fastq(.gz).
```bash
$linear read.fa genome.fa
``` 
Please add argument <b>'x'</b> if computing more than one reads and genomes.
```bash
$linear *fa x *fa
``` 
Use -h for more details of options
```bash
$linear -h
```

## Format of files
Output of Linear is an extended SAM/BAM format based on the standard formats for sequence alignment.
Definitions of each fields are extended to adapt to the alignment-free results.
Theses changes include:
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

- The 6th column of cigar is changed in the extended SAM/BAM.
Alignment cigar denotes each base, while alignment-free cigar denotes virtual alignment between 2 points.
Specially, alignment-free cigar in Linear is in the format of 'MG', where 'M' is allowed to be 'X' and '=' of standard cigar while 'G' is allowed to be 'I' and 'D' of standard cigar.
An example of virtual alignment of '200=80D' between 2 vertexes is shown in the following figure.
The green region is the estimated region, where the real alignment is supposed to be.

<img src="images/cigar_apx_map.png" alt="drawing" width="300"/>

- The 10th column of SEQ is inferred according to the reference and the 4,6th column rather than segment of read.
For cigar operation '=', the corresponding base from the reference rather than the read is inserted into the SEQ.
Thus the operation of '=' in alignment-free result doesn't necessarily mean the read is identical to the reference at the level of base pairs.
This is different from the SEQ for alignment.
The following table lists bases in SEQ from read or reference for each cigar operation, where
1/0 are Yes/No and 0.5 is conditional Yes

|From|=|X|M|I|D|S|H|
|--|--|--|--|--|--|--|--|
|read|0|0.5|0|1|0|1|0|
|ref|1|0|1|0|0|0|0|

- The 12th column, in which the definition of 'SA:Z' is changed because of the change of 6th and 10th columns.
Other tags are identical to the standard.
Standard definition of the tag can be found at [SAM/BAM format](https://samtools.github.io/hts-specs/SAMv1.pdf) and [Optional tags](https://samtools.github.io/hts-specs/SAMtags.pdf)

### Approximate mapping file (.apf)
This is a nonstandard format to provide readable overview of the alignment-free result.
The file can be enabled/disabled with the option '<b>-ot</b>'.
The .apf file contains the header and records
The following are the definition of the header and record.
### header
|col |filed|Description|Type|
|--|--|--|--|
|1|@|sign of header|{'@'}|
|2| QNAME|Query template NAME|string|
|3| QLEN|Query template LENGTH|int|
|4| QSTR|Query template mapped START| int |
|5| QEND|Query template mapped END| int |
|5| QSTRD|Query template mapped main STRAND|{'+','-'}|
|6| RNAME | Reference sequence NAME|String| 
|6| RLEN | Reference sequence LENGTH|int| 
|7| RSTR | Reference sequence mapped START|int| 
|8| REND | Reference sequence mapped END|int| 
### record
|col |filed|Description|Type|
|--|--|--|--|
|1|\||sign to start record|string|
|2|QSTR|Query region start|int|
|3|RSTR|Reference region start|int|
|4|DY|Distance of current 3th col to last|int|
|5|DX|Distance of current 4th col to last|int|
|6|RSTRD|region strand|{'+','-'}|

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

## Adaption to existing pipelines
The result of alignment is can be called by existing pipelines like other aligners.
Besides, the alignment-free results (with the alignment completely disbaled) can be called by existing alignment based pipelines, such as the SVs caller as well.

### Adaption to SVs callers
The compatibility of the alignment-free results to the SVs caller [PBSV](https://github.com/PacificBiosciences/pbsv) has been tested with PacBio long reads.

### Adaption to seqeunce graphical tools (IGV)
The compatibility of the alignment-free results to the IGV has been tested.
Alignment-free .bam  can be called directly.
