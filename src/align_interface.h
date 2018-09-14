#ifndef LINEAR_HEADER_ALIGNER_H
#define LINEAR_HEADER_ALIGNER_H

#include <seqan/align_parallel.h>

using namespace seqan;

int align (StringSet<String<Dna5> > & genomes,
           String<Dna5> & read, 
           String<Dna5> & comrevRead,
           String<uint64_t> & cord)
{
    typedef Align<String<Dna5>, ArrayGaps> TAlign; 
    typedef Row<TAlign>::Type TRow; 
    //using TExecPolicy = typename TestFixture::TExecPolicy;
    
    //TExecPolicy execPolicy;
    TAlign aligner; 
    resize(rows(aligner), 2); 
    //TRow & row1 = row(aligner, 0);
    //TRow & row2 = row(aligner, 1);
    double t, t1 = 0, t2 = 0;
    Infix<String<Dna5> >::Type infix1;  
    Infix<String<Dna5> >::Type infix2;  
    //String<Dna5> comrevRead;
    //_compltRvseStr(read, comrevRead);
    for (unsigned k = 1; k < length(cord) - 1; k++)
    {
        t = sysTime();
        uint64_t genomeId = _getSA_i1(_DefaultCord.getCordX(cord[k]));
        uint64_t genomeStart = _getSA_i2(_DefaultCord.getCordX(cord[k]));
        uint64_t strand = _DefaultCord.getCordStrand (cord[k]);
        uint64_t readStart = _DefaultCord.getCordY(cord[k]);
        //std::cout << "[]::align score " << genomeStart << " " << readStart << " " << genomeId << "\n";
        if (strand)
        {
            infix2 = infix(comrevRead, readStart, std::min(readStart + 192, length(read)));  
        }
        else
        {
            infix2 = infix(read, readStart, std::min(readStart + 192, length(read)));  
        }
        //infix1 = infix(genomes[genomeId], genomeStart, genomeStart + length(infix2));   
        infix1 = infix(genomes[genomeId], genomeStart, genomeStart + 192);   
        //assignSource (row(aligner, 0), infix1);  
        //assignSource (row(aligner, 1), infix2);  
        //clearClipping (row1);
        //clearClipping (row2);
        /*
        setClippedBeginPosition(row1, genomeStart);
        setClippedEndPosition(row1, genomeStart + 192);
        setClippedBeginPosition(row2, readStart);
        setClippedEndPosition(row2, readStart + 192);
        */
        t1 += sysTime() - t;
        t = sysTime();
        //int score = globalAlignment(aligner,  Score<int, Simple>(1, 0, 0), AlignConfig<false, false, false, false>(), -90, 90);
        int score = globalAlignmentScore(infix1,  infix2, Score<int, Simple>(1, 0, 0), AlignConfig<false, false, false, false>(), -90, 90);
        t2 += sysTime() - t;
        
        //std::cout << "[]::align score " << strand << " " << score << "\n" ;//<< aligner << "\n";
    }
    //std::cout << t1 << " " << t2 << " " << t1/t2 << std::endl;
    return 0;
}

int align(StringSet<String<Dna5> >& genomes, 
          StringSet<String<Dna5> > & reads, 
          StringSet<String<uint64_t> > & cords)
{
    double time = sysTime();

#pragma omp parallel
{
    #pragma omp for
    for (unsigned k = 0; k < length(reads); k += 1)
    {
        
        if (!empty(cords))
        {
            std::cerr << "read " << k << "\r";
            //align(genomes, reads[k], cords[k]);
        }
    }
}
    
    std::cerr << "[]::align time " << sysTime() - time << "\n";
    return 0;
}
#endif 
