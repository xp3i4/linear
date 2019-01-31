#ifndef LINEAR_HEADER_ALIGNER_H
#define LINEAR_HEADER_ALIGNER_H

//#include <seqan/align_parallel.h>

using namespace seqan;
/**
 * seqan::view coordinate and source coordinate transformation NOTES
 * setClippedBeginPosition(row1, c1)
 * @c1 := current_view_coordinate + clippedBeginPosition(row1)
 */
//TODO:: change score type
int const s1 = 3; //match
int const s2 = 0; //mismatch
int const s3 = -1; //gap
Score<int, Simple> score1_(s1, s2, s3);
float s_score_density_thd = 2; //if < the value alignment of cords will be dropped
int thd_align_score = 350 /*depends on score_scheme*/;

typedef Align<String<Dna5>, ArrayGaps> TAlign;
typedef Row<TAlign>::Type TRow; 
typedef Iterator<TRow>::Type TRowIterator;

uint64_t emptyCord = ~0;
/*
 * Structure of recods of start and end coordinates with aligner of gaps
 * Each record of c_pairs has two rows in the a_rows;
 * Gaps are defined as pairs of its start and end cords.
 * NOTE::structure for align_cords() function, different from gaps defined in 'gaps.h';
 */
class GapRecords
{
    typedef Row<Align<String<Dna5>, ArrayGaps> >::Type TRow;
    typedef std::pair<uint64_t, uint64_t> TCPair;
    typedef std::pair<TRow, TRow> TRPair;
    typedef String<std::pair<uint64_t, uint64_t> > TCPairs;
    typedef String<std::pair<TRow, TRow> > TRPairs;

public:
    TCPairs c_pairs;  //pair of start and end coordinate of gaps 
    TRPairs r_pairs;  //pair of alingments of blocks contains the start and end of gaps 
    String<int> bam_segs_id;
    uint64_t dx, dy;  //shift from start and end of gaps to the start and end of of the alignment

    TCPair & get_c_pair(int i) {return c_pairs[i];}  
    TRPair & get_r1_pair(int i){return r_pairs[i * 2];} //get rows of first(start) cord of i_th gaps
    TRPair & get_r2_pair(int i){return r_pairs[i * 2 + 1];} //get rows of second(end) cord of i_th gaps
    uint64_t getJointHeadCord(int i) 
    {
        return _DefaultCord.shift(c_pairs[i].first, -dx, -dy);
    }
    uint64_t getJointTailCord(int i) 
    {
        return _DefaultCord.shift(c_pairs[i].second, -dx, -dy);
    }
    int clear_(){
        clear(c_pairs); 
        clear(r_pairs);
        return 0;
    }
};

class ClipRecords: public GapRecords
{
public:
    uint64_t headCord; 
    uint64_t tailCord;
};
void setClippedPositions(TRow & row1, TRow & row2, int beginPos, int endPos)
{
    setClippedBeginPosition(row1, beginPos);
    setClippedBeginPosition(row2, beginPos);
    setClippedEndPosition(row1, endPos);
    setClippedEndPosition(row2, endPos);
}
/**
 *  return score of two chars
 */
//static float ln[10] = {1,2,}
inline int getScore_(char r1_char,
                     char r2_char,
                     int k,
                     int & x
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
inline int getScore_(TRowIterator & it1,
                     TRowIterator & it2,
                     int k,
                     int & x
)
{
    int s_match = 8 * k;
    int s_ins = -3 * k;
    int s_del = -3 * k;
    int s_mismatch = 0 * k;
    if (*it1 == *it2)
    {
        x += s_match;
    }
    else
    {
        if (isGap(it1))
        {
            x += s_ins;
        }
        else
        {
            if (isGap(it2))
                x += s_del;
            else
                x += s_mismatch;
        }
    }
    return x;
}
inline int getScore2_(TRowIterator & it1,
                     TRowIterator & it2,
                     int & len,
                     int & type,
                     int k,
                     int & x
)
{
    int s_match = 2 * k;
    int s_ins = -1 * k;
    int s_del = -1 * k;
    int s_mismatch = 0 * k;
    if (*it1 == *it2)
    {
        x += s_match * len;
    }
    else
    {
        if (isGap(it1))
        {
            x += s_ins;
        }
        else
        {
            if (isGap(it2))
                x += s_del;
            else
                x += s_mismatch;
        }
    }
    return x;
}
int inline setBamRecord (BamAlignmentRecord & bam_record,
                    Row<Align<String<Dna5>, ArrayGaps> >::Type & row1,
                    Row<Align<String<Dna5>, ArrayGaps> >::Type & row2, 
                    int g_id,
                    int g_beginPos
                   )
{
    align2cigar(bam_record.cigar, row1, row2);
    bam_record.rID = g_id;
    bam_record.beginPos = g_beginPos; 
}
/**
 * debug utility
 */
void printAlign_(Align<String<Dna5>, ArrayGaps> & aligner, int row_i, int row_j)
{
	std::cout << "[]::printAlignment \n";
    TRow & row1 = row(aligner, row_i);
    TRow & row2 = row(aligner, row_j);
    CharString line0, line1, line2, line3, line4;
    int line_len = 50;
    int sourceP = clippedBeginPosition(row1), sourceP2 = clippedBeginPosition(row2);
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
    for (int i = clippedBeginPosition(row1); i < clippedEndPosition(row1); i++)
    {
        int count = i % line_len;
        if (i == viewP)
        {
            CharString tmpc;
            int tmp = sourceP;
            while (tmp > 0)
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
        if (count == line_len - 1 || i + 1== length(row1))
        {
            std::cout << "     " << line0 << "\n     " << line1 << "\n     " << line3 << "\n     " << line2 << "\n     " << line4 << "\n";
            for (int j = 0; j < length(line0); j++)
            {
                line0[j] = ' ';
            }
			for (int j = 0; j < length(line1); j++)
            {
                line1[j] = ' ';
                line2[j] = ' ';
                line3[j] = ' ';
            }
            for (int j = 0; j < length(line4); j++)
            {
                line4[j] = ' ';
            }
       }
    }
}
void printAlignment(Align<String<Dna5>, ArrayGaps> & aligner)
{
	printAlign_(aligner, 0, 1);
}

inline int align_block (Row<Align<String<Dna5>, ArrayGaps> >::Type & row1,
                        Row<Align<String<Dna5>, ArrayGaps> >::Type & row2,
                        String<Dna5> & genome,
                        String<Dna5> & read,
                        String<Dna5> & comrevRead,
                        uint64_t strand,
                        uint64_t genomeStart,
                        uint64_t genomeEnd,
                        uint64_t readStart,
                        uint64_t readEnd,
                        int band
                       )
{
    //std::cout << "align len " << readStart << " " << readEnd << "\n";
    Infix<String<Dna5> >::Type infix1;  
    Infix<String<Dna5> >::Type infix2;  
    if (strand)
    {
        infix2 = infix(comrevRead, readStart, std::min(readEnd, length(read)));  
    }
    else
    {
        infix2 = infix(read, readStart, std::min(readEnd, length(read)));  
    }
    infix1 = infix(genome, genomeStart, std::min(genomeEnd, length(genome)));   
    assignSource (row1, infix1);  
    assignSource (row2, infix2); 
    int score = globalAlignment(row1, row2, Score<int, Simple> (s1, s2, s3), AlignConfig<true, true, true, true>(), -band, band);
    return score;
}

int cord2row_(Row<Align<String<Dna5>, ArrayGaps> >::Type & row1,
                Row<Align<String<Dna5>, ArrayGaps> >::Type & row2,
                String<Dna5> & genome,
                String<Dna5> & read, 
                String<Dna5> & comrevRead,
                uint64_t & cord,
                int block_size = window_size,
                int band = window_size / 2
               )
{
    uint64_t genomeStart = _getSA_i2(_DefaultCord.getCordX(cord));
    uint64_t genomeEnd = genomeStart + block_size;
    uint64_t readStart = _DefaultCord.getCordY(cord);
    uint64_t readEnd = readStart + block_size;
    uint64_t strand = _DefaultCord.getCordStrand (cord);   
    Infix<String<Dna5> >::Type infix1;  
    Infix<String<Dna5> >::Type infix2;  
    if (strand)
    {
        infix2 = infix(comrevRead, readStart, std::min(readEnd, length(read)));  
    }
    else
    {
        infix2 = infix(read, readStart, std::min(readEnd, length(read)));  
    }
    infix1 = infix(genome, genomeStart, std::min(genomeEnd, length(genome)));       
    assignSource (row1, infix1);  
    assignSource (row2, infix2); 
    return 0;
}

int align_cord (Row<Align<String<Dna5>, ArrayGaps> >::Type & row1,
                Row<Align<String<Dna5>, ArrayGaps> >::Type & row2,
				String<Dna5> & genome,
                String<Dna5> & read, 
                String<Dna5> & comrevRead,
                uint64_t & cord,
                int block_size = window_size,
                int band = window_size / 2
               )
{
    cord2row_ (row1, row2, genome, read, comrevRead, cord, block_size, band);
    int score = globalAlignment (row1, row2, Score<int, Simple> (s1, s2, s3), AlignConfig<true, true, true, true>(), -band, band);
    return score;
}

/**
 * Clip head or tails \in [0, l) of the aligner within the lxl window.
 * This function is to clip the well aligned region in the aligner.
 */
int clip_head_(Row<Align<String<Dna5>, ArrayGaps> >::Type & row1,
               Row<Align<String<Dna5>, ArrayGaps> >::Type & row2,
			   int g_end //view coordinates
			  )
{
    int window = 6;   //sliding window size
    int thd_clip = 5; //5 matches in the window
    int x = 0;
    int maxx = 0, maxxp = 0;
    int clip_start = clippedBeginPosition(row1);
    int flag = 0;
    int shift_head_len = -1;
    TRowIterator it1 = begin(row1);
    TRowIterator it2 = begin(row2);
    TRowIterator it1_2 = it1, it2_2 = it2; 
    if (clip_start > g_end - clip_start)
    {
    	return 1;
    }
    for (int i = 0; i < window; i++)
    {
   		if (*it1 == *it2)
   		{
   			++x;
   		}
   		it1++;
   		it2++;
    }
    for (int k = 0; k < g_end; k++)
    {
        if (x >= thd_clip)
    	{
    		setClippedBeginPosition(row1, k);
    		setClippedBeginPosition(row2, k);
            return 0;
    	}
    	if (*it1 == *it2)
    	{
    		++x;
    	}
    	if (*it1_2 == *it2_2)
    	{
            if (!flag)
            {
                shift_head_len = k;     //get fist match as len of free gap at the head
                flag = 1;
            }
    		--x;
    	}
    	if (maxx < x)
    	{
    		maxx = x;
    		maxxp = k;
    	}
    	++it1;
    	++it2;
    	++it1_2;
    	++it2_2;
    }
	setClippedBeginPosition(row1, maxxp);
	setClippedBeginPosition(row2, maxxp);
    return 0;
}

int clip_tail_(Row<Align<String<Dna5>, ArrayGaps> >::Type & row1,
               Row<Align<String<Dna5>, ArrayGaps> >::Type & row2,
			   int g_start
			  )
{
    TRowIterator it1, it2, it1_2, it2_2;
    int window = 6;
    int thd_clip = 5;
    int x = 0;
    int maxx = 0, maxxp = 0;
    int clip_start = clippedBeginPosition(row1);
    int clip_end = clippedEndPosition(row1);
    int flag = 0;
    int shift_tail_len = -1;
    it1 = end(row1) - 1;
   	it2 = end(row2) - 1;
   	it1_2 = it1;
   	it2_2 = it2;
    if (clip_end < g_start - clip_start)
    {
    	return 1;
    }
    for (int i = 0; i < window; i++)
    {
    	if (*(it1) == (*it2))
   		{
   			++x;
   		}
   		it1--;
   		it2--;
    }
    for (int k = clip_end; k > g_start; k--)
    {
    	if (x >= thd_clip)
    	{
    		setClippedEndPosition(row1, k);
    		setClippedEndPosition(row2, k);
    		return 0;
    	}
    	if (*it1 == *it2)
    	{
    		++x;
    	}
    	if (*it1_2 == *it2_2)
    	{
    		--x;
    	}
    	if (maxx < x)
    	{
    		maxx = x;
    		maxxp = k;
    	}
    	--it1;
    	--it2;
    	--it1_2;
    	--it2_2;
    }
    setClippedEndPosition(row1, maxxp);
    setClippedEndPosition(row2, maxxp);
    return 0;
}

int drop_align_(Row<Align<String<Dna5>, ArrayGaps> >::Type & row1,
                 Row<Align<String<Dna5>, ArrayGaps> >::Type & row2,
                 int score_align
                 //int g_end //view coordinates
                )
{
    int b1 = beginPosition(row1), b2 = beginPosition(row2);
    int e1 = endPosition(row1), e2 = endPosition(row2);
    float score_density = (float)score_align / std::max(e2 - b2, e1 - b1);
    if (score_density < s_score_density_thd) 
    {
        return 1; 
    }
    return 0;
}
/**
 * Merge any size of two blocks which starts from cord1 and cord2.
 * @parm rows are supposed to contain the alignment of the two blocks.
 * @parm rows are compared and merged from the corresponding clippedBeginPosition() to clippedEndPosition() 
 */
int merge_align__(Row<Align<String<Dna5>,ArrayGaps> >::Type & row11,
				 Row<Align<String<Dna5>,ArrayGaps> >::Type & row12,
				 Row<Align<String<Dna5>,ArrayGaps> >::Type & row21,
				 Row<Align<String<Dna5>,ArrayGaps> >::Type & row22,
				 uint64_t & cord1,  
				 uint64_t & cord2,
                 std::pair<int, int> & clips
				)
{
    if (cord1 == emptyCord)
    {
        return 0;
    }
    int thd_cord_overlap = std::max(endPosition(row11), endPosition(row12));
    if (!_DefaultCord.isCordsOverlap(cord1, cord2, thd_cord_overlap) ||
        _DefaultCord.getCordStrand(cord1 ^ cord2))  //sv: gap or reverse
    {
        std::cout << "[]::merge_align_ " 
                  << thd_cord_overlap << " "
                  << get_cord_x(cord1) << " "
                  << get_cord_x(cord2) << " "
                  << "\n";
        return 1;
    }
    int bit = 20, bit2 = 40;
    uint64_t start11 = _getSA_i2(_DefaultCord.getCordX(cord1));
    uint64_t start21 = _getSA_i2(_DefaultCord.getCordX(cord2));
    uint64_t start12 = _DefaultCord.getCordY(cord1);
    uint64_t start22 = _DefaultCord.getCordY(cord2);
    int64_t mask = (1ULL << bit) - 1;
	int64_t delta1 = start21 - start11;
	int64_t delta2 = start22 - start12;
	String<int64_t> align1, align2; 
    TRowIterator it1, it2;
    if (endPosition(row11) < beginPosition(row21) + delta1 ||
        endPosition(row12) < beginPosition(row22) + delta2 ||
        endPosition(row11) > endPosition(row21) + delta1 ||
        endPosition(row12) > endPosition(row22) + delta2 
       )
    {
        std::cout << "[]xxx " 
                  << endPosition(row11) << " " 
                  << endPosition(row12) << " " 
                  << endPosition(row21) << " " 
                  << beginPosition(row21) << " "
                  << delta1 << " "
                  << endPosition(row22) << " "
                  << beginPosition(row22) << " "
                  << delta2 << " "
                  << "\n";
        return 2;
    }
    //ouput source coordinates of align to buffer 
    int64_t src1 = beginPosition(row21) + delta1;  //cord1_x
    int64_t intersect_view_Begin = toViewPosition(row11, src1);
    int64_t src2 = toSourcePosition(row12, intersect_view_Begin); 
    int64_t intersect_view_End = clippedEndPosition(row11) - clippedBeginPosition(row11);
    it1 = begin(row11) + intersect_view_Begin;
    it2 = begin(row12) + intersect_view_Begin;
	for (int64_t i = intersect_view_Begin; i < intersect_view_End; i++)
	{
		if (*it1 == *it2)
		{
			appendValue (align1, ((i + clippedBeginPosition(row11))<< bit2) + (src1 << bit) + src2);
		}
        if (!isGap(it1))
        {
            src1++;
        }
        if (!isGap(it2))
        {
            src2++;
        }
		it1++; 
        it2++;
	}

    src1 = beginPosition(row21); //cord1_x
    src2 = beginPosition(row22);
    intersect_view_Begin = 0;

    //intersect_view_End = toViewPosition(row21, endPosition(row11) - delta1);
    it1 = begin(row21) + intersect_view_Begin;
    it2 = begin(row22) + intersect_view_Begin;
    for (int64_t i = intersect_view_Begin; i < intersect_view_End; i++)
    {
        if (*it1 == *it2)
        {
            appendValue (align2, ((i + clippedBeginPosition(row21))<< bit2) 
                                 + ((src1 + delta1) << bit) + (src2 + delta2));
        }
        if (!isGap(it1))
        {
            src1++;
        }
        if (!isGap(it2))
        {
            src2++;
        }
        it1++; 
        it2++;
    }
    if (length(align1) == 0 || length(align2) == 0)
    {
        return 3;
    }
    
    int thd_merge_x = 2, thd_merge_y = 2;
    int flag = 0, start_j = 0;
    int64_t x1 = (align1[0] >> bit) & mask;
    int64_t y1 = align1[0] & mask;flag = 0;
    int sum = 0;
	for (int i = 0; i < length(align1) - 1; i++)	
	{
		int64_t x1_next = (align1[i + 1] >> bit) & mask;
		int64_t y1_next = align1[i + 1] & mask;
        flag = 0;
		for (int j = start_j; j < length(align2); j++)	
	    {
            sum++;
			int64_t x2 = align2[j] >> bit & mask;
			int64_t y2 = align2[j] & mask;
            if (!flag)
            {
                if (std::abs(x1_next - x2) < thd_merge_x) 
                {
                    start_j = j;
                    flag = 1;
                }
            }
			if (std::abs(x1 - x2) < thd_merge_x && std::abs(y1 - y2) < thd_merge_y)
			{
                int clip1 = (align1[i] >> bit2 & mask);
                int clip2 = (align2[j] >> bit2 & mask);
                clips.first = clip1; clips.second = clip2;
                std::cout << "[]::merge_align_ " << sum << "\n";
				return 0;
			}
			else if (x2 - x1 > thd_merge_x|| y2 - y1 > thd_merge_x)
			{
				break;
			}
		}
        x1 = x1_next;
        y1 = y1_next;
	}
    return 4;
}
int merge_align_(Row<Align<String<Dna5>,ArrayGaps> >::Type & row11,
                 Row<Align<String<Dna5>,ArrayGaps> >::Type & row12,
                 Row<Align<String<Dna5>,ArrayGaps> >::Type & row21,
                 Row<Align<String<Dna5>,ArrayGaps> >::Type & row22,
                 uint64_t & cord1,  
                 uint64_t & cord2
                )
{
    std::pair<int, int> clips;
    
    int flag = merge_align__(row11, row12, row21, row22, cord1, cord2, clips);
    if (!flag)
    {
        int clip1 = clips.first;
        int clip2 = clips.second;
        setClippedBeginPosition(row21, clip2);
        setClippedBeginPosition(row22, clip2);
        setClippedEndPosition(row11, clip1);
        setClippedEndPosition(row12, clip1);
    }
    return flag;
}

inline void insertGaps(GapRecords & gaps,
                  uint64_t cord1,    //start coordinate of gap
                  uint64_t cord2,    //end coordinate of gap
                  int bam_segs_id,
                  int thd_merge_gap,
                  int dx_,
                  int dy_
                 )
{
    if (empty(gaps.c_pairs))
    {
        appendValue(gaps.c_pairs, std::pair<uint64_t, uint64_t>(cord1, cord2));
        appendValue(gaps.bam_segs_id, bam_segs_id);

    }
    else 
    {
        if (_DefaultCord.isCordsOverlap(back(gaps.c_pairs).second, cord1, thd_merge_gap) && 
            back(gaps).bam_segs_id == bam_segs_id)
        {
            int i = length(gaps.c_pairs) - 1;
            gaps.get_c_pair(i).second = cord2;
        }
        else
        {
            appendValue(gaps.c_pairs, std::pair<uint64_t, uint64_t>(cord1, cord2));
            appendValue(gaps.bam_segs_id, bam_segs_id);
        }
    }
    gaps.dx = dx_;
    gaps.dy = dy_;
}

inline void insertGaps(GapRecords & gaps,
                  uint64_t cord1,    //start coordinate of gap
                  uint64_t cord2, // end coordinate
                  Row<Align<String<Dna5>, ArrayGaps> >::Type & row11,
                  Row<Align<String<Dna5>, ArrayGaps> >::Type & row12,
                  Row<Align<String<Dna5>, ArrayGaps> >::Type & row21,
                  Row<Align<String<Dna5>, ArrayGaps> >::Type & row22,
                  int bam_segs_id,
                  int thd_merge_gap,
                  int dx_,
                  int dy_ 
                 )
{
    if (empty(gaps.c_pairs))
    {
        appendValue(gaps.c_pairs, std::pair<uint64_t, uint64_t>(cord1, cord2));
        appendValue(gaps.r_pairs, std::pair<TRow, TRow>(row11, row12));
        appendValue(gaps.r_pairs, std::pair<TRow, TRow>(row21, row22));
        appendValue(gaps.bam_segs_id, bam_segs_id);
    }
    else 
    {
        if (_DefaultCord.isCordsOverlap(back(gaps.c_pairs).second, cord1, thd_merge_gap) && 
            back(gaps).bam_segs_id == bam_segs_id)
        {
            int i = length(gaps.c_pairs) - 1;
            gaps.get_c_pair(i).second = cord2;
            gaps.get_r2_pair(i).first = row21;
            gaps.get_r2_pair(i).second = row22;
        }
        else
        {
            appendValue(gaps.c_pairs, std::pair<uint64_t, uint64_t>(cord1, cord2));
            appendValue(gaps.r_pairs, std::pair<TRow, TRow>(row11, row12));
            appendValue(gaps.r_pairs, std::pair<TRow, TRow>(row21, row22));
            appendValue(gaps.bam_segs_id, bam_segs_id);
        }
    }
    gaps.dx = dx_;
    gaps.dy = dy_;
}
inline void _nextView2Src (Iterator<Row<Align<String<Dna5>, ArrayGaps>>::Type>::Type & it1,
                    Iterator<Row<Align<String<Dna5>, ArrayGaps>>::Type>::Type & it2, 
                    int64_t & src1, 
                    int64_t & src2)
{
    if (!isGap(it1))
    {
        src1++;
    }
    if (!isGap(it2))
    {
        src2++;
    }
}
/**
 *  The size of the block to be clipped are supposed to > window_size,
 *  since the function clips the middle part of the block excluding 
 *  the head and tail part. Too small block cann't be well clipped. 
 */
int clip_segs_(Row<Align<String<Dna5>,ArrayGaps> >::Type & row1,
               Row<Align<String<Dna5>,ArrayGaps> >::Type & row2,
               Row<Align<String<Dna5>,ArrayGaps> >::Type & row_jh1,
               Row<Align<String<Dna5>,ArrayGaps> >::Type & row_jh2,
               Row<Align<String<Dna5>,ArrayGaps> >::Type & row_jt1,
               Row<Align<String<Dna5>,ArrayGaps> >::Type & row_jt2,
               String<std::pair<int, int> > & clip_records,
               uint64_t gap_start_cord,
               uint64_t cord_jh,
               uint64_t cord_jt,
               Score<int> & score_scheme,
               int thd_clip_score,
               int thd_reject_score,
               int thd_accept_score,
               int thd_accept_density,
               int thd_clip_mini_region = 20,
               int window = 30, //supposed to < 50
               int delta = 3
              )
{
    int x = 0;
    int thd_1 = 1;
    int flag = 0;
    int64_t src1 = beginPosition(row1);
    int64_t src2 = beginPosition(row2);
    TRowIterator it1 = begin(row1);
    TRowIterator it2 = begin(row2);
    TRowIterator it1_2 = it1, it2_2 = it2; 
    String<int> buffer;  //score at the intervals by delta
    String<int> buffer_src1; 
    String<int> buffer_src2; //alignment source coordiante
    int clipped_end = clippedEndPosition(row1); 
    int clipped_begin = clippedBeginPosition(row1);
    int clipped_len = clipped_end - clipped_begin;
    if (clipped_len < window * 3)
    {
        return 1;
    }
    int buf_len = 0;
    int count = delta;
    resize (buffer, clipped_len);
    resize (buffer_src1, clipped_len);
    resize (buffer_src2, clipped_len);
    int mth_len = 0, ins_len = 0, del_len = 0, mis_len = 0;
    for (int k = 0; k < clipped_len; k++)
    {
        getScore_(it1, it2, 1, x);
        _nextView2Src(it1, it2, src1, src2);
        if (count++ == delta)
        {
            buffer[buf_len] = x; //+ thd_1 * (std::abs(src1 - src1_2 - src2 + src2_2));
            buffer_src1[buf_len] = src1;
            buffer_src2[buf_len] = src2;
            //std::cout << "buffer " << buffer_src2[buf_len] << " " << buf_len << "\n";
            buf_len++;
            count = 1;
        }
        ++it1;
        ++it2;
    }
//TODO::change the score so it will not clip long ins or dels.
    resize (buffer, buf_len);
    int clip_beginPos = 0; 
    int clip_endPos = 0;
    int last_prebp = 0;
    int prebp = 0; //pre-breakpoint counter of buffer
    int bp = prebp; //bp * delta = view coordinate
    int bp_sc1 = beginPosition(row1), bp_sc2 = beginPosition(row2); //breakpoint of source 
    int last_accept_prebp_sc1, last_accept_prebp_sc2;
    int last_accept_bp_view, last_accept_prebp_view;
    int prebp_sc1 = bp_sc1, prebp_sc2 = bp_sc2;
    int d_w = window / delta;
    int d_m = thd_clip_mini_region / delta ;
//TODO!!::clip the fist and last segment
    int last_ = length(buffer) - std::max(d_w, d_m);
    int ct_clips = 0;  //count number of clip records
    std::pair<int, int> clip_pair;
    for (int i = d_w + 1; i < last_ ; i++)
    {
        int d_score_left = buffer[i] - buffer[i - d_w];
        int d_score_right = buffer[i + d_w] - buffer[i];
        int d_score = buffer[i + d_w] - (buffer[i] << 1) + buffer[i - d_w];
        //std::cout << "d_score " << i << " " << d_score  << " " << buffer[i + d_w] - buffer[i] << " " << buffer[i] - buffer[i - d_w] << "\n";
        if ((std::abs(d_score) > thd_clip_score && 
             std::min(d_score_left, d_score_right) < thd_reject_score)|| i == last_ - 1) 
        {
            bp = i;
            int max_d_score = std::abs(d_score);
            int d_j = std::min(i + d_m, (int)length(buffer) - d_w);
            for (int j = i + 1; j < d_j; j++)
            {
                //std::cout << "d_score " << j << " " << d_score  << " " << buffer[j + d_w] - buffer[j] << " " << buffer[j] - buffer[j - d_w]  << " " << bp << " " << max_d_score<< "\n";
                d_score = buffer[j + d_w] - (buffer[j] << 1) + buffer[j - d_w]; 
                if (std::abs(d_score) > max_d_score)
                {
                    bp = j;
                    max_d_score = std::abs(d_score);
                } 
            } //get the maximal clip_score in region specified thd_clip_mini_region 
              //when the first d_score satisifed appear
            bp_sc1 = buffer_src1[bp];
            bp_sc2 = buffer_src2[bp];


//TODO::score to accecpt the segs, need to make score more reasonable
            if ((buffer[bp] - buffer[prebp]) / (bp - prebp + 1) > thd_accept_density)
            {
                std::cout << "[]::clip_segs_ " 
                          << prebp * delta << " "
                          << bp * delta << " "
                          << bp_sc1 << " "
                          << bp_sc2 << " "
                          << (buffer[bp] - buffer[prebp])/ (bp - prebp + 1) << " "
                          << toViewPosition(row1, bp_sc1) << " "
                          << toViewPosition(row2, bp_sc2) << " "
                          << get_cord_x(gap_start_cord) << " "
                          << get_cord_y(gap_start_cord) << "\n"; 
                clip_beginPos = prebp * delta + clipped_begin;
                clip_endPos = bp * delta + clipped_begin;
                clip_pair.first = clip_beginPos;
                clip_pair.second = clip_endPos;
                if (ct_clips == 0)
                {
                    uint64_t headCord = _DefaultCord.shift(gap_start_cord, prebp_sc1, prebp_sc2);
                    setClippedPositions(row1, row2, clip_beginPos, clip_endPos);
                    flag = merge_align_(row_jh1, row_jh2, row1, row2, cord_jh, headCord);
                }
                appendValue(clip_records, clip_pair);
                last_accept_prebp_sc1 = prebp_sc1;
                last_accept_prebp_sc2 = prebp_sc2;
                last_accept_bp_view = bp;
                last_accept_prebp_view = prebp; 
                ++ct_clips;
            }
            prebp_sc1 = bp_sc1;
            prebp_sc2 = bp_sc2; 
            prebp = bp;
            i += d_j - i;
        }
    }
    clip_beginPos = last_accept_prebp_view * delta + clipped_begin;
    clip_endPos = last_accept_bp_view * delta + clipped_begin; 
    clip_pair.first = clip_beginPos;
    clip_pair.second = clip_endPos;
    uint64_t tailCord = gap_start_cord;
    setClippedPositions(row1, row2, clip_beginPos, clip_endPos);
    flag = merge_align_(row1, row2, row_jt1, row_jt2, tailCord, cord_jt); 
    return 0;
}
int merge_segs_(String<BamAlignmentRecord> & bam_records,
               StringSet<String<std::pair<int, int> > > & clip_records,
               TAlign & aligner,
               String<uint64_t> & gap_start_cords
              )
{
    int count_ht = 0;
    int new_bam_size = 0;
    int bam_it = length(bam_records);
    for (int j = 0; j < length(clip_records); j++)
    {
        new_bam_size += length(clip_records[j]);
    }
    //resize (bam_records, new_bam_size);
    BamAlignmentRecord bam_record;
    for (int j = 0; j < length(clip_records); j++)
    {
        for (int k = 0; k < length(clip_records[j]); k++)
        {
            setClippedPositions(row(aligner, 2 * j), row (aligner, 2 * j + 1), clip_records[j][k].first, clip_records[j][k].second);
            uint64_t g_id = _getSA_i1(_DefaultCord.getCordX(gap_start_cords[j]));
            uint64_t g_beginPos = get_cord_x (gap_start_cords[j]) + beginPosition(row (aligner, 2 * j));
       
            setBamRecord(bam_record, row(aligner, 2 * j), row(aligner, 2 * j + 1), g_id, g_beginPos);
            appendValue(bam_records, bam_record);
            std::cout << "[]::merge_segs " 
                     << clip_records[j][k].first << " " 
                     << clip_records[j][k].second << " "
                     << "\n";
            for (int m = 0; m < length(bam_record.cigar[m].operation); m++)
            {
                std::cout << bam_record.cigar[m].count << bam_record.cigar[m].operation;
            }

            ++bam_it;
        }
    }
    return 0;
}
int align_gaps (String<BamAlignmentRecord> & bam_records,
                GapRecords & gaps,
                StringSet<String<Dna5> >& genomes,
                String<Dna5> & read, 
                String<Dna5> & comrevRead,
                Score<int> & score_scheme,
                int thd_clip_score,
                int thd_reject_score,
                int thd_accept_score,
                int thd_accept_density
               )
{
    if (empty(gaps.c_pairs))
    {
        return 0;
    }
    std::cout << "[]align_gaps\n ";
    Align<String<Dna5>, ArrayGaps> aligner;
    StringSet<String<std::pair<int, int> > > clip_records;
    String<uint64_t> gap_start_cords;
    String<uint64_t> gap_end_cords;
    resize(rows(aligner), 4);
    resize(clip_records, 2);
    resize(gap_start_cords, 2);
    resize(gap_end_cords, 2);
    int flag = 0;
    for (int i = 0; i < length(gaps.c_pairs); i++)
    {
        gap_start_cords[0] = gaps.get_c_pair(i).first;
        gap_end_cords[0] = gaps.get_c_pair(i).second;
        cmpRevCord (gap_start_cords[0], gap_end_cords[0], gap_start_cords[1], gap_end_cords[1], length(read));
        for (int j = 0; j < length(clip_records); j++)
        {
            clear (clip_records[j]);
            uint64_t gap_start_cord = gap_start_cords[j];
            uint64_t gap_end_cord = gap_end_cords[j];
            if (_DefaultCord.getCordStrand(gap_start_cord ^ gap_end_cord))
            {
                return 1;
            }
            int block_size = std::max(get_cord_x(gap_end_cord - gap_start_cord),
                                      get_cord_y(gap_end_cord - gap_start_cord));
            int g_id = _getSA_i1 (_DefaultCord.getCordX(gap_start_cord));
            int head_end = block_size >> 2;
            int tail_start = block_size - (block_size >> 2);
            int band = block_size >> 1;

            align_cord (row(aligner, 2 * j), row(aligner, 2 * j + 1), genomes[g_id], read, comrevRead, gap_start_cord, block_size,band);
            clip_segs_(row(aligner, 2 * j), 
                       row(aligner, 2 * j + 1), 
                       gaps.get_r1_pair(i).first,
                       gaps.get_r1_pair(i).second,
                       gaps.get_r2_pair(i).first,
                       gaps.get_r2_pair(i).second,
                       clip_records[j], 
                       gap_start_cord, 
                       gaps.getJointHeadCord(i),
                       gaps.getJointTailCord(i),
                       score_scheme, 
                       thd_clip_score,
                       thd_reject_score,
                       thd_accept_score,
                       thd_accept_density
                      );
        }
        merge_segs_(bam_records, 
                clip_records, 
                aligner,
                gap_start_cords
             //   gaps.get_r1_pair(i).first,
             //   gaps.get_r1_pair(i).second,
             //   gaps.get_r2_pair(i).first,
             //   gaps.get_r2_pair(i).second,
             //   gap_start_cord, 
             //   gaps.getJointHeadCord(i),
             //   gaps.getJointTailCord(i),
                );
    }
    return 0;
}

/*
set_cigar_soft_clip
                    //int n = length(read) * strand - _nStrand(strand) * (_DefaultCord.getCordY(cords[i]) + endPosition(row(aligner, ri + 1)));
                    //appendValue(back(bam_records).cigar, CigarElement<>('S', n));

                resize(bam_records, length(bam_records) + 1);
                //int n = length(read) * strand - _nStrand(strand) * (_DefaultCord.getCordY(cords[i]) + beginPosition(row(aligner, ri + 1)));
                //appendValue(back(bam_records).cigar, CigarElement<>('S', n));
                back(bam_records).flag = (back(bam_records).flag & (~16)) | (_DefaultCord.getCordStrand(cords[i]) << 4); 
*/
/*
 * debug utility
 */
void printGaps(String<std::pair<uint64_t, uint64_t> > & gaps)
{
    for (int i = 0; i < length(gaps); i++)
    {
        std::cout << "[]::gaps " 
                  << _getSA_i1(_DefaultCord.getCordX(gaps[i].first)) << " " 
                  << _getSA_i2(_DefaultCord.getCordX(gaps[i].first)) << " " 
                  << _DefaultCord.getCordY(gaps[i].first) << " " 
                  << _getSA_i1(_DefaultCord.getCordX(gaps[i].second)) << " " 
                  << _getSA_i2(_DefaultCord.getCordX(gaps[i].second)) << " " 
                  << _DefaultCord.getCordY(gaps[i].second) << "\n"; 
    }
}
/*
 * Align cords and output cigar string.
 * Each cord will be clipped if necessary. 
 */
int align_cords (StringSet<String<Dna5> >& genomes,
                 String<Dna5> & read, 
                 String<Dna5> & comrevRead,
                 String<uint64_t> & cords,
                 String<BamAlignmentRecord> & bam_records,
                 int block_size = window_size,
                 int band = window_size / 2
                ) 
{
    Align<String<Dna5>, ArrayGaps> aligner;
    GapRecords gaps;
    BamAlignmentRecord emptyBamRecord;
    String<uint64_t> cords_buffer;
    int head_end = block_size >> 2, tail_start = block_size - (block_size >> 2);
    int ri = 0, ri_pre = 2; //cliped segment and row id
    int g_id = -1, g_beginPos = 0, strand = 0, flag = 0;
    int thd_merge_gap = block_size; // two adjacent gaps will be merged to one if < it
    int flag2 = 0;
    uint64_t preCord = emptyCord; //last well aligend and accepted cord
    resize(cords_buffer, 2);
    resize(rows(aligner), 4); 
    double t1, t2 = 0, t3 = sysTime();
    if (length(cords) < 2) // cords is empty
    {
        return 0;
    }
    for (int i = 1; i < (int)length(cords); i++)
    {
        //>debug section begin xxxxxxx
        if (i == 10)
            cords[i] += 50;
        if (i > 44 && i < 53)
            cords[i] += 100;
        //<debug section end
        g_id = _getSA_i1(_DefaultCord.getCordX(cords[i]));
        g_beginPos = _getSA_i2(_DefaultCord.getCordX(cords[i]));
        strand = _DefaultCord.getCordStrand(cords[i]);
        flag = 0;
        int score_align = 
            align_cord (row(aligner, ri), 
                        row(aligner, ri + 1), 
                        genomes[g_id], 
                        read, 
                        comrevRead, 
                        cords[i]);
        flag |= clip_head_ (row(aligner, ri), row(aligner, ri + 1), head_end);
        flag |= clip_tail_ (row(aligner, ri), row(aligner, ri + 1), tail_start);
        flag |= drop_align_(row(aligner, ri), row(aligner, ri + 1), score_align);
        std::cout << "[]::align_cords_1 " << i << " " << flag << " " << get_cord_y(cords[i]) << " " << get_cord_x(cords[i]) <<"\n";
        if (flag)
        {
            continue; //alignment quality check, drop poorly aligned 
        }

        if (_DefaultCord.isBlockEnd(cords[i - 1]))
        {
            preCord = emptyCord;
        }
        flag = merge_align_(
                          row(aligner, ri_pre), 
                          row(aligner, ri_pre + 1),
                          row(aligner, ri),
                          row(aligner, ri + 1),
                          preCord, 
                          cords[i]
                         );
        if (flag) //merge failed 
        {
            int dx = block_size >> 1; //shift the cord to the gap cord;
            int dy = block_size >> 1;
            insertGaps(gaps, 
                       _DefaultCord.shift(preCord, dx, dy),
                       _DefaultCord.shift(cords[i], dx, dy),
                       row(aligner, ri_pre), 
                       row(aligner, ri_pre + 1),
                       row(aligner, ri),
                       row(aligner, ri + 1),
                       length(bam_records) - 1,
                       thd_merge_gap,
                       dx,
                       dy);
            //appendValue (bam_records, emptyBamRecord);
            flag2 = 1;
            std::cout << "[]::align_cords cord_i " 
                      << get_cord_x(gaps.getJointTailCord(length(gaps.c_pairs) - 1)) << " "
                      << get_cord_x(cords[i]) << " \n" 
                      << "\n";
        }
        /*
        else
        {
            if (!flag2) // when the rows_pre are not gap_tailJoint cord.
            {
                setBamRecord(back(bam_records), 
                         row(aligner, ri_pre), 
                         row(aligner, ri_pre + 1),
                         g_id,
                         g_beginPos
                        );
            }
            flag2 = 0;
        }
        */
        preCord = cords[i];
        //std::cout << "[]::align_cords_2 " << i << " " << flag << " " << preCord << " " << cords[i] << "\n";
        std::swap (ri, ri_pre); //swap the current and pre row id in the aligner.
    }
    printGaps(gaps.c_pairs);
    Score<int> score_scheme;
    int thd_clip_score = 80;
    int thd_reject_score = 130;
    int thd_accept_score = 140;
    int thd_accept_density = 16;
    align_gaps(bam_records, gaps, genomes, read, comrevRead, score_scheme, thd_clip_score, thd_reject_score, thd_accept_score, thd_accept_density); 
    return 0;
}
/*
void clip_cigar_(String<CigarElement<> > & cigars)
{
    for (unsigned i = 0; i < length(cigars); i++)
    {

    }
}
*/ 
/**
 * Clip breakpoint \in [w, l - w),w = 30, within the lxl window
 * Direction: if  > 0  -----------mmmmmmmmm,  if < 0 mmmmmmmmmmm--------; 
 * where 'm' represents well aligned part
 */
int clip_window_(Align<String<Dna5>,ArrayGaps> & aligner, 
				 int g_start,
				 int g_end,
			     uint64_t & clip_ref, 
				 uint64_t & clip_read, 
				 int direction
				)
{
	typedef Align<String<Dna5>, ArrayGaps> TAlign;
	typedef Row<TAlign>::Type TRow; 

	double t = sysTime();
	TRow & row1 = row(aligner, 0);
    TRow & row2 = row(aligner, 1);
    int window = 30;  // w = window
    int x = 0;
    String<int> buffer;
//WARNING & TODO::toViewPosition is extreamly time inefficient
    for (int i = toViewPosition(row1, g_start); i < toViewPosition(row1, g_start + window); i++)
    {
        x = getScore_ (row1[i], row2[i], 1, x);
    }
    int delta = 3;
    for (int k = delta + g_start; k < g_end - window  + 1; k += delta)
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
    }
    t = sysTime() - t;
    double t2 = sysTime();
    int max = 0;
    int max_sp_ref = 0; // max source position
    int max_sp_read = 0;
    direction = direction / std::abs(direction);
    for (int k = window / delta ; k < (int)length(buffer); k++)
    {
        int d_ = (buffer[k] - buffer[k - window / delta]) * direction;
        if (max < d_)
        {
            max = d_;
            max_sp_ref = k * delta;
            max_sp_read = toSourcePosition(row2, toViewPosition(row1, max_sp_ref));
        }
    }
    clip_ref = (max < 20)?-1:max_sp_ref;
    clip_read = (max < 20)?-1:max_sp_read;
    return 0;
}
/**
 * Clip break point of the alignment of genome and read within the lxl window 
 * direction: clip direction  > 0  -----------mmmmmmmmm,  < 0 mmmmmmmmmmm--------; 
 * where 'm' is match
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
	                         int direction,
	                         int score_flag = 1
	                        )
{ 
    typedef Row<Align<String<Dna5>,ArrayGaps> >::Type TRow;
    Align<String<Dna5>,ArrayGaps> aligner;
    resize(rows(aligner), 2); 
    TRow row1 = row(aligner, 0);
    TRow row2 = row(aligner, 1);
    for (int i = genomeStart; i < genomeEnd; i++)
    {
        std::cout << *(begin(genome) + i);
    }
    int score = align_block  (row1,
                              row2,
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
    if (score < thd_align_score / (int) window_size * g_range && score_flag)
    {
        return -1;
    }
    uint64_t clip_ref = 0, clip_read = 0;
  	clip_window_ (aligner, 0, g_range, clip_ref, clip_read, direction);
    uint64_t returnCord = _DefaultCord.createCord(_createSANode(genomeId, genomeStart + clip_ref), 
                                                  readStart + clip_read, 
                                                  strand);
    return returnCord;
}

#endif 
