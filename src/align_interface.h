#ifndef LINEAR_HEADER_ALIGNER_H
#define LINEAR_HEADER_ALIGNER_H

//#include <seqan/align_parallel.h>

using namespace seqan;


int const s1 = 3; //match
int const s2 = 0; //mismatch
int const s3 = -1; //gap

int thd_align_score = 350 /*depends on score_scheme*/;

/**
 * debug utility
 */
void printAlignment(Align<String<Dna5>, ArrayGaps> & aligner)
{
    typedef Align<String<Dna5>, ArrayGaps> TAlign;
    typedef Row<TAlign>::Type TRow;  
    TRow & row1 = row(aligner, 0);
    TRow & row2 = row(aligner, 1);
    CharString line0, line1, line2, line3, line4;
    int line_len = 50;
    int sourceP = 0, sourceP2 = 0;
    int viewP = toSourcePosition(row1, sourceP);
    int viewP2 = toSourcePosition(row2, sourceP);
    int len0 = 10;
    resize(line0, line_len * 2);
    resize(line1, line_len);
    resize(line2, line_len);
    resize(line3, line_len);
    resize(line4, line_len * 2);
    for (int j = 0; j < length(line0); j++)
    {
        line0[j] = ' ';
    } 
    for (int j = 0; j < length(line0); j++)
    {
        line4[j] = ' ';
    }
    for (int i = 0; i < length(row1); i++)
    {
        int count = i % line_len;
        if (i == viewP)
        {
            CharString tmpc;
            int tmp = sourceP;
            while ( tmp > 0)
            {
               appendValue(tmpc, tmp - tmp / 10 * 10 + '0');
               tmp /= 10;
            }
            for (int j = 0; j < length(tmpc); j++)
            {
                line0[count + length(tmpc) - 1 - j] = tmpc[j];
            //    std::cerr << line0[count + length(tmpc) - 1 - j];
            }
            //line0[count] = '0';
            sourceP += len0;
            viewP = toViewPosition(row1, sourceP - 1);
        }
        if (i == viewP2)
        {
            CharString tmpc;
            int tmp = sourceP2;
            while ( tmp > 0)
            {
               appendValue(tmpc, tmp - tmp / 10 * 10 + '0');
               tmp /= 10;
            }
            for (int j = 0; j < length(tmpc); j++)
            {
                line4[count + length(tmpc) - 1 - j] = tmpc[j];
            }
            sourceP2 += len0;
            viewP2 = toViewPosition(row2, sourceP2 - 1);
        }
        line1[count] = row1[i];
        line2[count] = row2[i];
        if (line1[count] == line2[count])
        {
            line3[count] = '|';
        }
        else
        {
            line3[count] = ' ';
        }
        if (count == line_len - 1)
        {
            std::cout << "     " << line0 << "\n     " << line1 << "\n     " << line3 << "\n     " << line2 << "\n     " << line4 << "\n";
            for (int j = 0; j < length(line0); j++)
            {
                line0[j] = ' ';
            }
            for (int j = 0; j < length(line4); j++)
            {
                line4[j] = ' ';
            }
        }
    }
}
/*
int align_merge_cord_band (StringSet<String<uint64_t> > & cords,
                           String<uint64_t> & bands,
                           int band_width,
                           int band_width_max
                          )
{
    int flag = 1;
    for (int i = 1; i < length(cords); i++)
    {
        int64_t tmp_upper_band;
        int64_t tmp_lower_band;
        if (flag)
        {
            int upper_band = _DefaultCord.getCordX(cords[i]) + band_width; 
            int lower_band = _DefaultCord.getCordY(cords[i]) + ;
            flag = 0;
        }
        else
        {
            int64_t x1 = _DefaultCord.getCordX(cords[i]); 
            int64_t y1 = _DefaultCord.getCordY(cords[i]);
            ///WARNING::need optimize
            if (x1 - upper_band - y1  + band_width> 0 )
            {
                tmp_upper_band = x1 - y1;
            }
            else if (x1 - lower_band - y1 +  band_width <= 0)
            {
                    tmp_lower_band = x1 - y1;
            }
            if (x1 - lower_band - y1  - band_width > 0 )
            {
                tmp_lower_band_band = x1 - y1;
            }
            else if (x1 - upper_band - y1 -  band_width <= 0)
            {
                    tmp_upper_band = x1 - y1;
            }
            if (tmp_upper_band - tmp_lower_band < band_width_max) 
            {
                upper_band = tmp_upper_band;
                lower_band = tmp_lower_band;
            }
            else // clip the bands to two groups of discontinuous bands otherwise the band width is too large.
            {
                uint64_t center_diagonal = lower_band + band_width;
                uint64_t x_start1 = 
                if (x_start - y_start  - center_diagonal > 0)
                {
                    
                }
                else
                {
                    x1 - y1 - 
                }
                appendValue(bands, start_cords);
                appendValue(bands, end_cords);
            }
        }
        if (_DefaultCord.isCordEnd(cords))
        {
            flag == 1;
        }
    }
}
*/
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
    double time = sysTime ();
    int score = globalAlignment(aligner, Score<int, Simple> (s1, s2, s3), AlignConfig<true, true, true, true>(), -band, band);
    double dt1 = sysTime() - time;
    int score2 = globalAlignmentScore(infix1, infix2, Score<int, Simple> (s1, s2, s3), AlignConfig<true, true, true, true>(), -band, band);
    double dt2 = sysTime() - time;
    std::cout << "[]::align_block_ " << dt1 / dt2 << "\n";
    printAlignment(aligner);
    std::cout << "[]::align_block " << " score " << score << " " << genomeStart << " " << genomeEnd << " " << readStart << " " << readEnd << " " << strand << " " << band << "\n" ;
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
    double t1 = sysTime();
    typedef Align<String<Dna5>,ArrayGaps> TAlign;
    typedef Row<TAlign>::Type TRow; 
    TAlign aligner; 
    resize(rows(aligner), 2); 
    double t2 = sysTime();
    for (int i = genomeStart; i < genomeEnd; i++)
    {
        std::cout << *(begin(genome) + i);
    }
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
    double dt = sysTime() - t1;
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
        //std::cout << "[]::clip_window " << d_ << " " << k << " " << buffer[k] << " " << buffer[k - window] << "\n";;
        
        if (max < d_)
        {
            max = d_;
            max_sp_ref = k * delta;
            max_sp_read = toSourcePosition(row2, toViewPosition(row1, max_sp_ref));
        }
        
    }
    clip_ref = (max < 20)?-1:max_sp_ref;
    clip_read = (max < 20)?-1:max_sp_read;
    //std::cout << "[]::clip_window::clip " << clip_ref << " " << clip_read << " " << genomeStart + clip_ref << " " << readStart + clip_read << "\n";
    
    uint64_t returnCord = _DefaultCord.createCord(_createSANode(genomeId, genomeStart + clip_ref), 
                                                  readStart + clip_read, 
                                                  strand);
    //std::cout << "[]::clip_window::max_score " << max << " " << int(score) / g_range << " " << _DefaultCord.getCordX(returnCord) << "\n";
    std::cout << "[]::clip_window time percent " << dt / (sysTime() - t1) << " " << (t2 - t1) / dt << "\n";
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
