## Linear

Mapping reads efficiently. 

## Usage

read length < 1MBases

## Prerequisites

| Platform                                              | Building tools            |
| ----------------------------------------------------- | ------------------------- |
| <img src="images/lx_icon.png" width="30"> Linux 4.9.0 | gcc  > 4.9; cmake > 3.0.0 |

## Install

```bash
  $ CMake .
  $ make linear 
```

## Built with

### Included in the source
- [SeqAn 2.0](<https://seqan.readthedocs.io/en/master/>)-Library for sequence analysis

- [googletest](<https://github.com/google/googletest>)-Unit test(optional)

### Pre-install 
- [zlib](<https://www.zlib.net/>)-(File compression I/O)
If not installed by default on Debian (Ubuntu, Mint,..) based distributions:
```bash
   $sudo apt install zlib1g-dev libbz2-dev
```

## File format
### Approximate mapping file (.amf)

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

  

## license

Linear is licensed under the 3-clause BSD license

## Contact



<img src="/home/cx/code/linear/linear/images/logo.svg" width="330">





