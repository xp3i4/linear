# qbdnb
 
Tools for mapping reads efficiently. 

## Usage

* Use as standalone app please clone the master branch to SeqAn2 master branch/apps

```bash
  $ CMake -DCMAKE_BUILD_TYPE=Release [directory to SeqAn]
  $ make pacMapper
```

* Use as library please include the mapper.h in the source.

## Interface 
Major classes are encapsulated as the `class Mapper` in mapper.h as part of the pipline. 

- Filter
  - Using `Mapper::hits()` to return the StringSet of best `hits` for each read 
  - Using `Mapper::getHitX(hit)` and `Mapper::getHitY(hit)` to return the coordinates of hit
- Band estimation
  - Using `Mapper::cords()` returns the StringSet of coordinates `cords` of slidings windows for each read.
  - Using `Mapper::getCordX(cord)` and `Mapper::getCordY(cord)` to return the coordinates of sliding window
- map
  - To add new module to the pipeline please provide interfaces in th `Mapper` class and call them in the `map()` function 
  or I will encapsulated them at last.
  
## Notice

The current version is to complete the interface and the pipeline. While the benchmark should be conducted later after all the things are accomplished.  







