# qbdnb

Prototype of tools for mapping smrt sequencing reads efficiently. 

Usa
To use it please note "include mapper.h" in the source file. 
and call map(genome, reads) in the pacmapper.h to return the regions
It takes StringSet\<Dna5\> as input, and return Pair<uint64_t, uint64_t> as the begin and end position in the genome 

# notice

The result is supposed to output

readsid, genomeid, length(reads), length(genome) begin, end, strand, score

1.For convenience it currently output begin and end .
please take the region: [begin, end] ± 15% * readlength in the genome for verification.

2.Revese complement strand dna is calculated during the process but the output section for them is not included yet. So please use forward strand data set for test.



