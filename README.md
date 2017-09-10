# qbdnb
 
Tools for mapping reads efficiently. 

## Usage

* Use as standalone app please clone the master branch to SeqAn2 master branch/apps

```bash
  $ CMake -DCMAKE_BUILD_TYPE=Release [directory to SeqAn]
  $ make pacMapper
```

* Use as library please #include the mapper.h in the source.

## Interface 
Classes and functions are encapsulated as the `class Mapper` in mapper.h so it can be easier for developing different modules. Other modules like the I/O will be add later. Please add interfaces for other module in the Mapper class so other modules can call them or I will encapsulated the them at last.

- Filter
  - Using `Mapper::hits()` to return the best `hits` for each read. 
  - Using `Mapper::getHitX(hit)` and `Mapper::getHitY(hit)` to return the coordinates of hit
- Band estimation
  - `Mapper::cords()` returns the coordinates `cords` of slidings windows for each read, Its return type is StringSet
  - `Mapper::getCordX(cord)` and `Mapper::getCordY(cord)` return the coordinates of sliding window

## Notice

The current version is to complete the interface and the pipeline. While the benchmark should be conducted later after all the development and tunning are accomplished.  







