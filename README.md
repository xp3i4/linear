#qbdnb

 
Prototype for mapping reads efficiently. 

Usage
--------

 Use standalone app please clone the master branch to SeqAn2 master branch/apps/

::

  $ CMake -DCMAKE_BUILD_TYPE=Release [directory to SeqAn]
  
  $ make pacMapper

Use classes and functions please #include the mapper.h in the source.
 Filter
Using mapper.hits() to return the best hits of each read. 
Using mapper.getHitX() and getHitY() to return the coordinates of hit
Band estimation
mapper.cords() returns the slidings windows of each read, Its return type is StringSet
mapper getCordX() and getCordY return the coordinates of sliding window

Notice
--------
Uniform interface is defined as the class Mapper in mapper.h to make it easier for developing and integerating different modules.  
Please add interfaces for other module in the Mapper class, so they can be called easilly. 

The repository is to provide the pipeline so other moduls can be add in or movd out. While benchmark should be conducted later.  







