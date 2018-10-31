#ifndef LINEAR_HEADER_ALIGNER_H
#define LINEAR_HEADER_ALIGNER_H

//#include <seqan/align_parallel.h>

using namespace seqan;


int const s1 = 3; //match
int const s2 = 0; //mismatch
int const s3 = -1; //gap

int thd_align_score = 350 /*depends on score_scheme*/;

int align_genome_cord (String<Dna5> & genome,
                       String<Dna5> & read, 
                       String<Dna5> & comrevRead,
                       String<uint64_t> & cord,
                       String<int> & score)
{
    Infix<String<Dna5> >::Type infix1;  
    Infix<String<Dna5> >::Type infix2;  
    for (unsigned k = 0; k < length(cord); k++)
    {
        uint64_t genomeStart = _getSA_i2(_DefaultCord.getCordX(cord[k]));
        uint64_t strand = _DefaultCord.getCordStrand (cord[k]);
        uint64_t readStart = _DefaultCord.getCordY(cord[k]);
        //std::cout << "[]::align score " << genomeStart << " " << readStart << " " << genomeId << "\n";
        if (strand)
        {
            infix2 = infix(comrevRead, readStart, std::min(readStart + window_size, length(read)));  
        }
        else
        {
            infix2 = infix(read, readStart, std::min(readStart + window_size, length(read)));  
        }
        infix1 = infix(genome, genomeStart, genomeStart + window_size);   

        score[k] = globalAlignmentScore(infix1, infix2, Score<int, Simple> (s1, s2, s3), AlignConfig<false, false, false, false>(), -90, 90);
        
    }
    return 0;
}

inline int align_block_(Align<String<Dna5>, ArrayGaps> & aligner,
                      String<Dna5> & genome,
                      String<Dna5> & read,
                      String<Dna5> & comrevRead,
                      uint64_t strand,
                      uint64_t genomeStart,
                      uint64_t genomeEnd,
                      uint64_t readStart,
                      uint64_t readEnd,
                      int band)
{

    Infix<String<Dna5> >::Type infix1;  
    Infix<String<Dna5> >::Type infix2;  

    if ( strand)
    {
        infix2 = infix(comrevRead, readStart, std::min(readEnd, length(read)));  
    }
    else
    {
        infix2 = infix(read, readStart, std::min(readEnd, length(read)));  
    }
    infix1 = infix(genome, genomeStart, genomeEnd);   
    assignSource (row(aligner, 0), infix1);  
    assignSource (row(aligner, 1), infix2); 
    int score = globalAlignment(aligner, Score<int, Simple> (s1, s2, s3), AlignConfig<true, true, true, true>(), -band, band);
    std::cout << "[]::align_block " << aligner << " score " << score << " " << genomeStart << " " << genomeEnd << " " << readStart << " " << readEnd << " " << strand << " " << band << "\n" ;
    return score;
}

/**
 * a short function to return score x of two chars
 */
//static float ln[10] = {1,2,}
inline int getScore_(char r1_char,
                     char r2_char,
                     int k,
                     //int & l,
                     int x
)
{
    int s_match = 8 * k;
    int s_gap = -3 * k;
    int s_mismatch = 0 * k;
    if (r1_char == r2_char)
    {
        x += s_match;
    }
    else
    {
        if (r1_char == '-')
        {
            x += s_gap;
        }
        else
        {
            x += s_mismatch;
        }
        
    }
    return x;
}

/**
 * Clip break point within the lxl window 
 * direction: clip direction  > 0  -----------mmmmmmmmm,  < 0 mmmmmmmmmmm--------; where 'm' is match
 */
inline uint64_t clip_window (String<Dna5> & genome,
                        String<Dna5> & read,
                        String<Dna5> & comrevRead,
                        uint64_t genomeId,
                        uint64_t genomeStart,
                        uint64_t genomeEnd,
                        uint64_t readStart,
                        uint64_t readEnd,
                        uint64_t strand,
                        uint64_t band,
                        int direction)
{ 
    typedef Align<String<Dna5>,ArrayGaps> TAlign;
    typedef Row<TAlign>::Type TRow; 
    
    TAlign aligner; 
    resize(rows(aligner), 2); 
    int score = align_block_ (aligner,
                              genome, 
                              read, 
                              comrevRead,
                              strand, 
                              genomeStart, 
                              genomeEnd,
                              readStart,
                              readEnd,
                              band
                            );
    int g_range = (int) genomeEnd - genomeStart;
    if (score < thd_align_score / (int) window_size * g_range)
    {
        return -1;
    }
    
    /**
     * clip breakpoint within the block
     */
    TRow & row1 = row(aligner, 0);
    TRow & row2 = row(aligner, 1);
    int window = 30;
    int x = 0;
    int count = 0;
    int count2 = 0;
    int s_match = 4;
    int clip_thd = 70;
    uint64_t clip_ref = 0;
    uint64_t clip_read = 0;
    String<int> buffer;
    x = 0;
    for (int i = toViewPosition(row1, 0); i < toViewPosition(row1, window); i++)
    {
        x = getScore_ (row1[i], row2[i], 1, x);
    }
    int delta = 3;
    for (int k = delta; k < g_range - window  + 1; k += delta)
    {
        for (int i = toViewPosition(row1, k - delta); i < toViewPosition(row1, k); i++)
        {
            x = getScore_ (row1[i], row2[i], -1, x);
        }
        for (int i = toViewPosition(row1, k + window - 1 - delta); i < toViewPosition(row1, k + window - 1); i++)
        {
            x = getScore_ (row1[i], row2[i], 1, x);
        }
        appendValue(buffer, x);
        //std::cout << "clip_block_right " << k << " " << x << "\n";
    }
    int max = 0;
    int max_sp_ref = 0; // max source position
    int max_sp_read = 0;
    direction = direction / std::abs(direction);
    for (int k = window / delta ; k < (int)length(buffer); k++)
    {
        int d_ = (buffer[k] - buffer[k - window / delta]) * direction;
        /**
        for (int i = toViewPosition(row1, k * delta - 10); i < toViewPosition(row1, k * delta + 10); i++) 
        {
            std::cout << row1[i];
        }
        std::cout << "\n";
        */
        std::cout << "[]::clip_window " << d_ << " " << k << " " << buffer[k] << " " << buffer[k - window] << "\n";;
        
        if (max < d_)
        {
            max = d_;
            max_sp_ref = k * delta;
            max_sp_read = toSourcePosition(row2, toViewPosition(row1, max_sp_ref));
        }
        
    }
    clip_ref = (max < 20)?-1:max_sp_ref;
    clip_read = (max < 20)?-1:max_sp_read;
    std::cout << "[]::clip_window::clip " << clip_ref << " " << clip_read << " " << genomeStart + clip_ref << " " << readStart + clip_read << "\n";
    
    uint64_t returnCord = _DefaultCord.createCord(_createSANode(genomeId, genomeStart + clip_ref), 
                                                  readStart + clip_read, 
                                                  strand);
    std::cout << "[]::clip_window::max_score " << max << " " << int(score) / g_range << " " << _DefaultCord.getCordX(returnCord) << "\n";
    return returnCord;
    //return clip_ref;
}

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
    std::string cigar;
    std::string mutations;
    for (unsigned k = 1; k < length(cord) - 1; k++)
    {
        t = sysTime();
        uint64_t genomeId = _getSA_i1(_DefaultCord.getCordX(cord[k]));
        uint64_t genomeStart = _getSA_i2(_DefaultCord.getCordX(cord[k]));
        uint64_t strand = _DefaultCord.getCordStrand (cord[k]);
        uint64_t readStart = _DefaultCord.getCordY(cord[k]);
        clear(cigar);
        clear(mutations);
        //std::cout << "[]::align score " << genomeStart << " " << readStart << " " << genomeId << "\n";
        if (strand)
        {
            infix2 = infix(comrevRead, readStart, std::min(readStart + window_size, length(read)));  
        }
        else
        {
            infix2 = infix(read, readStart, std::min(readStart + window_size, length(read)));  
        }
        //infix1 = infix(genomes[genomeId], genomeStart, genomeStart + length(infix2));   
        infix1 = infix(genomes[genomeId], genomeStart, genomeStart + window_size);   
        assignSource (row(aligner, 0), infix1);  
        assignSource (row(aligner, 1), infix2);  
        //clearClipping (row1);
        //clearClipping (row2);
        /*
        setClippedBeginPosition(row1, genomeStart);
        setClippedEndPosition(row1, genomeStart + window_size);
        setClippedBeginPosition(row2, readStart);
        setClippedEndPosition(row2, readStart + window_size);
        */
        t1 += sysTime() - t;
        t = sysTime();
        int score = globalAlignment(aligner, Score<int, Simple> (s1, s2, s3), AlignConfig<false, false, false, false>(), -90, 90);
        //int score = globalAlignmentScore(infix1,  infix2, score_scheme, AlignConfig<false, false, false, false>(), -90, 90);
        t2 += sysTime() - t;
        
        //std::cout << "[]::align score " << genomeStart << " " << readStart << " " << strand << " " << score << "\n" ;//<< aligner << "\n";
        //std::cout << "[]::align strand " << strand <<  " " << genomeStart << " score " << score << "\n" << aligner << "\n";
        align2cigar_(aligner, cigar, mutations);
        std::cout << "[]::align cigar " << cigar << " " << mutations << "\n";
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
