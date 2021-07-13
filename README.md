# Linear
## alignment-free methods 
Mapping reads efficiently. 

## Prerequisites

| Building tools  |   Version          |
| ------------------- | ------------------------- |
| Linux | 4.9.0|
| gcc|>4.9|
| cmake|>3.0.0|


### External libraries in the source

- [SeqAn 2.0](<https://seqan.readthedocs.io/en/master/>)-Library for sequence analysis

- [googletest](<https://github.com/google/googletest>)-Unit test

### Pre-install required
- [zlib](<https://www.zlib.net/>)-(File compression I/O).
  To install on Ubuntu based distributions:
```bash
   $sudo apt install zlib1g-dev libbz2-dev
```

## Build
Create a new build directory. In the build directory
```bash
  $ CMake [path to git cloned source] 
  $ make linear 
```

## Usage
Support.fa(.gz), .fastq(.gz) for input.
```bash
$linear read.fa genome.fa
``` 
Please add 'x' when mapping more than one reads and genomes.
```bash
$linear *fa x *fa
``` 
Use -h for more details regarding options
```bash
$linear -h
```

## File format of results
### SAM/BAM
Standard format for alignment and map.
The definition of SAM/BAM of alignment of Linear is identical to the standard.
The definition of SAM/BAM of map of Linear is changed as the following table.


|col |filed|Description|Status|
|--|--|--|--|
|   1  | QNAME | Query template NAME                       | Yes       |           
|   2  | FLAG  | bitwise FLAG                              | Yes       | 
|   3  | RNAME | Reference sequence NAME                   | Yes       | 
|   4  | POS   | 1-based leftmost mapping POSition         | Yes       | 
|   5  | MAPQ  | MAPping Quality                           | Yes       | 
|   6  | CIGAR | CIGAR string                              | Changed   | 
|   7  | RNEXT | Reference name of the mate/next read      | Yes       |
|   8  | PNEXT | Position of the mate/next read            | Yes       |
|   9  | TLEN  | observed Template LENgth                  | Yes       | 
|   10 | SEQ   | segment SEQuence                          | Changed   |
|   11 | QUAL  | ASCII of Phred-scaled base QUALity+33     | Yes       |
|   12 | TAG   | Optional tags                             | Changed   |

- 6th column of cigar

<img src="images/cigar_apx_map.png" alt="drawing" width="300"/>
- 10th column of SEQ is inferred according to the reference and the 4,6th column rather than segment of read.
For cigar operation '=', the corresponding base from the reference rather than the read is inserted into the SEQ.
Thus the operation of '=' in result of mapping doesn't necessarily mean the read is identical to the reference at the level of base pairs.
This is different from the SEQ for alignment.
The change of definitation is to make the SEQ of mapping compatible to existing tools
For operations of 'M', 'X', 'I', and 'S' the corresponding bases in the read are inserted.
This is identical to the SEQ for alignment.
- 12th column, in which the definition of 'SA:Z' is changed because of the change of 6th and 10th columns.
Other tags are identical to the standard.
Standard definition of the tag can be found at [SAM/BAM format](https://samtools.github.io/hts-specs/SAMv1.pdf) and [Optional tags](https://samtools.github.io/hts-specs/SAMtags.pdf)

### Approximate mapping file (.apf)
This is a file

```
@> S1_14253 10049 45 10049 + chr1 44948 9964 19212 
| 45 9964 0 0 1 
| 215 10119 170 155 2 
| 352 10240 137 121 3 
| 446 10327 94 87 4 
...
| 2400 12112 144 128 17 
| 2544 12256 144 144 18 
| 2688 12384 144 128 19 
| 2832 12528 144 144 20 
| 2826 12515 -6 -13 21 
|**+++++++++++ 6786 12607 3960 92 22 
| 3118 12771 -3668 164 23 
| 3196 12845 78 74 24 
...
| 6656 16080 144 128 48 
| 6800 16208 144 128 49 
| 6944 16336 144 128 50 
|**+++++++++++ 1931 16507 -5013 171 51 
|**+++++++++++ 2126 16690 195 183 52 
|**+++++++++++ 2309 16857 183 167 53 
|**+++++++++++ 2462 16999 153 142 54 
|**+++++++++++ 2592 17120 130 121 55 
|**+++++++++++ 2674 17189 82 69 56 
...
| 9680 18839 0 7 68 
| 9870 19020 190 181 69 
```
[sam]()-sam file
[gvf]()

## Adapation to SVs callers
Some of existing tools uses character of space rather than the tab in the SAM/BAM header. 
Though it violates the standard definition, they are widly 
Enable the output of .bam for pbsv with the option of <b>-o 8</b>.
## Adaption to seqeunce graphical tools (IGV)

## license
BSD License 2.0


## Contact









