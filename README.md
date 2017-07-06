# qbdnb

this reprository provoides a mapper prototype for smrt sequencing reads. 
the structures will be polished soon.

to use it please include the header files 
and call map(genome, reads) in the pacmapper.h to return the regions
It takes StringSet\<Dna5\> as input, and return Pair<uint64_t, uint64_t> as the begin and end position in the genome 

please take the region: [begin, end] ± 15% * readlength in the genome for verification

