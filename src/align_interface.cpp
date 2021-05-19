#include <utility> 
#include "base.h"
#include "cords.h" 
#include "align_util.h"
//#include "f_io.h"
#include "align_bands.h"
#include "align_interface.h"
//#include "new_funcs.h"
//TODO seqand::setclippedpositoin retrieve source postion that's not efficient
//TODO make holder to rows, so set clip postion can be iteratored.
using namespace seqan;
/**
 * Notation regarding the Seqan::(clip function) declaration that might cause confusion
 * clippedBeginPosition() returns "unclipped view pos" 
 * setClippedBeginPosition(gaps, c1); closed [c1,)
 * insertGap(gap, pos) pos is the "clipped! view pos"
 * setclippedEndPosition(gaps, c1)  open[, c1)
 * @c1 := current_view_coordinate + clippedBeginPosition(row1)
 */

unsigned _default_block_size_ = 96; // make sure the value == window_size in pmpfinder.h
/**
 * debug util
 */
void printAlign_(TRow & row1, TRow & row2)
{
    std::cout << "[]::printAlignment \n";
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
    printAlign_(row(aligner, 0), row(aligner, 1));
}
void printCigars(String<CigarElement< > > &cigar, std::string header)
{
    std::cout << header;
    for (int i = 0; i < length(cigar); i++)
    {
        std::cout << cigar[i].count <<cigar[i].operation;
    }
    std::cout << "\n";
}

class GapRecordHolder{
    GapRecords & holder;
    int it;
public:
    GapRecordHolder(GapRecords & gaps);

    int getBamSegIdHead();
    int getBamSegIdTail();
    int atEnd();
    GapRecordHolder & next();
    GapRecordHolder & operator [] (int);
    GapRecordHolder & operator ++();
    typename GapRecords::TCPair & getCords();
};
GapRecordHolder::GapRecordHolder(GapRecords & gaps):
    holder(gaps), it(0){}
inline int GapRecordHolder::atEnd()
{
    return it >= length(holder.c_pairs);
}
inline int GapRecordHolder::getBamSegIdTail()
{
    return holder.getBamSegIdTail(it);
}
inline int GapRecordHolder::getBamSegIdHead()
{
    return holder.getBamSegIdHead(it);
}
inline GapRecordHolder & GapRecordHolder::next()
{
    it++; 
    return *this;
}
inline GapRecordHolder & GapRecordHolder::operator [] (int i)
{
    it = i;
    return *this;
}
inline GapRecordHolder & GapRecordHolder::operator ++()
{
    return next();
}
inline typename GapRecords::TCPair & 
GapRecordHolder::getCords()
{
    return holder.get_c_pair(it);
}

//TODO:: change score type
int const s1 = 3; //match
int const s2 = -2; //mismatch
int const s3 = -1; //gap extend
int const s4 = -1; //gap open
float s_score_density_thd = 2; //if < the value alignment of cords will be dropped
float s_score_window_thd = 0.75;

uint64_t CORD_NULL = _DefaultCord.makeBlockEndVal(~0);
uint64_t emptyCord = CORD_NULL - 1;

AlignGapParms _gap_parm;
Score<int, Simple> _default_scheme_ (s1, s2, s3, s4);
typename GapRecords::TCPair & GapRecords::get_c_pair(int i) {return c_pairs[i];}  
typename GapRecords::TRPair & GapRecords::get_r1_pair(int i) 
{
    if (i < 0)
    {
        return r_pairs[length(r_pairs) - 2];
    }
    return r_pairs[i * 2];
}
typename GapRecords::TRPair & GapRecords::get_r2_pair(int i)
{
    if (i < 0)
    {
        return r_pairs[length(r_pairs) - 1];
    }
    return r_pairs[i * 2 + 1];
} 
int GapRecords::getBamSegIdHead (int i){return bam_segs_id[i];}
int GapRecords::getBamSegIdTail (int i){return bam_segs_id[i] + 1;}
int GapRecords::get_clip_flag (int i)
{
    if (empty(clip_flags))
    {
        return flag_clip_unset;
    }
    if (i < 0)
    {
        return back(clip_flags);
    }
    return clip_flags[i];
}
int GapRecords::set_clip_flag(int clip_flag, int pos)
{
    if (empty(clip_flags))
        return 1;
    int p = pos;
    if (pos < 0 || p >= length(clip_flags))
    {
        p = length(clip_flags) - 1;
    }
    int value = clip_flag ;
    clip_flags[p] = value;
    return 0;
}
uint64_t GapRecords::getJointHeadCord(int i) 
{
    if (empty (c_pairs))
    {
        return CORD_NULL;
    }
    if (i == -1 || i >= length(c_pairs))
    {
        return _DefaultCord.shift(back(c_pairs).first, -dx, -dy);
    }
    return _DefaultCord.shift(c_pairs[i].first, -dx, -dy);
}
uint64_t GapRecords::getJointTailCord(int i) 
{
    if (empty (c_pairs))
    {
        return CORD_NULL;
    }
    if (i == -1 || i >= length(c_pairs))
    {
        return _DefaultCord.shift(back(c_pairs).second, -dx, -dy);
    }
    return _DefaultCord.shift(c_pairs[i].second, -dx, -dy);
}
int GapRecords::clear_(){
    clear(c_pairs); 
    clear(r_pairs);
    return 0;
}
AlignGapParms::AlignGapParms ():
        thd_clip_score(80),   //density score: clip when density is < 0.8
        thd_reject_score(130),      //calculate density within the window
        thd_accept_score(140),  //the mini distanse between two clips.
        thd_min_interval(20),
        thd_accept_density(4.5),
        cbrlht_thd_src_d_bps(200),
        thd_clip_view_len(1000)
{
    thd_clip_scheme = Score<int> (1, -1, -1, -1);
}
void setClippedBeginPositions(TRow & row1, TRow & row2, int beginPos)
{
    setClippedBeginPosition(row1, beginPos);
    setClippedBeginPosition(row2, beginPos);
}           
void setClippedEndPositions(TRow & row1, TRow & row2, int endPos)
{
    setClippedEndPosition(row1, endPos);
    setClippedEndPosition(row2, endPos);
}
void setClippedPositions(TRow & row1, TRow & row2, int beginPos, int endPos)
{
    setClippedBeginPositions(row1, row2, beginPos);
    setClippedEndPositions(row1, row2, endPos);
}
/**
 *  return score of two chars
 */
//static float ln[10] = {1,2,}
 int getScore_(char r1_char,
               char r2_char,
               int k,
               int & x
               )
{
    int s_match = 8 * k;
    int s_gap = -3 * k;
    int s_mismatch = -3 * k;
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
 int getScore_(TRowIterator & it1,
                     TRowIterator & it2,
                     int k,
                     int & x
)
{
    int s_match = 8;
    int s_ins = -3;
    int s_del = -3;
    int s_mismatch = 0;
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

 int cumulate_match_sc__ (int l, int s)
{
    int j = 5;
    int maxk = 5;
    return std::min((l / j + 1), maxk) * s;
}
 int cumulate_gap_sc__ (int l, int s)
{
    int j = 2;
    int maxk = 5;
    return std::min((l / j + 1), maxk) * s;
}
/*
 * @flag:=[1]bit|[31]type
 * bit:= 0 continuse, 1 seg breakpoint
 * type:= match 0, mismatch 2, ins 4, del 8, gap open 16
 */
 int getScore2_(TRowIterator & it1,
                     TRowIterator & it2,
                     Score<int, Simple> score_scheme,
                     int & l1, //match
                     int & l2, //mismatch
                     int & l3, //ins
                     int & l4, //del
                     int & l5, //gap_open
                     int & x,
                     int & sege_flag,
                     int & type_flag
)
{
    int s1 = score_scheme.data_match; 
    int s2 = score_scheme.data_mismatch;
    int s3 = score_scheme.data_gap_extend;
    int s4 = score_scheme.data_gap_open;
    int thd_continuous_gap = 2;
    if (*it1 == *it2)
    {
        x += cumulate_match_sc__(l1, s1);
        l1++;
        l2 = 0;
        l5 = 0;
        if (type_flag != 1)
        {
            sege_flag = 1;
        }
        else
        {
            sege_flag = 0;
        }
        type_flag = 1;
    }
    else
    {
        if (isGap(it1))
        {
            x += cumulate_gap_sc__(l3, s4);
            l1 = l2 = l4 = l5 = 0;
            l3++;
            if (type_flag != 2)
            {
                type_flag = 2; 
                sege_flag = 1;
            }
            else
            {
                sege_flag = 0;
            }
        }
        else
        {
            if (isGap(it2))
            {
                x += cumulate_gap_sc__(l4, s4);
                l1 = l2 = l5 = l3 = 0;
                l4++;
                if (type_flag != 4)
                {
                    type_flag = 4; 
                    sege_flag = 1;
                }
                else
                {
                    sege_flag = 0;
                }                
            }
            else
            {
                x += s2;
                l1 = l5 = 0;
                l2++;
                if (type_flag != 8)
                {
                    type_flag = 8; 
                    sege_flag = 1;
                }
                else
                {
                    sege_flag = 0;
                }    
            }
        }
    }
    return x;
}
 int getScore2_(TRowIterator & it1,
                     TRowIterator & it2,
                     Score<int, Simple> score_scheme,
                     String<int> & ls,
                     int & x,
                     int & sege_flag,
                     int & type_flag)
{
    if (length(ls) < 5)
    {
        resize (ls, 5);
        for (int i = 0; i < length(ls); i++)
        {
            ls[i] = 0;
        }
    }
    return getScore2_(it1, it2, score_scheme, ls[0], ls[1], ls[2], ls[3], ls[4], x, sege_flag, type_flag);
}

/**
 * Replace operator '=' of row
 * Otherwise sequece of row1 depends on sequence of row2. see seqan::Holder
 */
int detach_copy_row(TRow & row1, TRow & row2)
{
    row1 = row2;
    detach(row1);
}

int align_block (Row<Align<String<Dna5>, ArrayGaps> >::Type & row1,
                 Row<Align<String<Dna5>, ArrayGaps> >::Type & row2,
                 String<Dna5> & genome,
                 String<Dna5> & read,
                 String<Dna5> & comrev_read,
                 uint64_t strand,
                 uint64_t genomeStart,
                 uint64_t genomeEnd,
                 uint64_t readStart,
                 uint64_t readEnd,
                 int band,
                 Score<int> scheme = _default_scheme_)
{
    //std::cout << "align len " << readStart << " " << readEnd << "\n";
    Infix<String<Dna5> >::Type infix1;  
    Infix<String<Dna5> >::Type infix2;  
    if (strand)
    {
        infix2 = infix(comrev_read, readStart, std::min(readEnd, length(read)));  
    }
    else
    {
        infix2 = infix(read, readStart, std::min(readEnd, length(read)));  
    }
    infix1 = infix(genome, genomeStart, std::min(genomeEnd, length(genome)));   
    assignSource (row1, infix1);  
    assignSource (row2, infix2); 
    int score = globalAlignment(row1, row2, scheme, AlignConfig<true, true, true, true>(), -band, band);
    return score;
}

int cord2row_(Row<Align<String<Dna5>, ArrayGaps> >::Type & row1,
              Row<Align<String<Dna5>, ArrayGaps> >::Type & row2,
              String<Dna5> & genome,
              String<Dna5> & read, 
              String<Dna5> & comrev_read,
              uint64_t cord_str,
              uint64_t cord_end)
{
    uint64_t genomeStart = get_cord_x(cord_str);
    uint64_t genomeEnd = get_cord_x(cord_end);
    uint64_t readStart = get_cord_y(cord_str);
    uint64_t readEnd = get_cord_y(cord_end);
    uint64_t strand = get_cord_strand (cord_str);   
    if (genomeStart > genomeEnd || readStart > readEnd)
    {
        return 1;
    }
    Infix<String<Dna5> >::Type infix1;  
    Infix<String<Dna5> >::Type infix2;  
    if (strand)
    {
        infix2 = infix(comrev_read, readStart, std::min(readEnd, length(read)));  
    }
    else
    {
        infix2 = infix(read, readStart, std::min(readEnd, length(read)));  
    }
    infix1 = infix(genome, genomeStart, std::min(genomeEnd, length(genome)));       
    //std::cout << "rowsize " << length(infix1) << " " << length(infix2) << "\n";
    assignSource (row1, infix1);  
    assignSource (row2, infix2); 
    return 0;
}
int align_cord (Row<Align<String<Dna5>, ArrayGaps> >::Type & row1,
                Row<Align<String<Dna5>, ArrayGaps> >::Type & row2,
                String<Dna5> & genome,
                String<Dna5> & read, 
                String<Dna5> & comrev_read,
                uint64_t cord_str,
                uint64_t cord_end,
                int band_lower,
                int band_upper,
                int local_flag = 1,
                Score<int> scheme = _default_scheme_)
{
    if (cord_str == cord_end || 
        get_cord_x(cord_str) == get_cord_x(cord_end) || 
        get_cord_y(cord_str) == get_cord_y(cord_end))
    {
        return 0;
    }
    cord2row_ (row1, row2, genome, read, comrev_read, cord_str, cord_end);
    int score = 0;
    if (!local_flag)
    {
        score = localAlignment (row1, row2, scheme, DynamicGaps());
    }
    else
    {
     //   double time = sysTime();
        score = globalAlignment (row1, row2, scheme, AlignConfig<true, true, true, true>(), -band_lower, band_upper);
    //    std::cout << "atime " << sysTime() - time << "\n";
    }
    return score;
}
int align_cord (Row<Align<String<Dna5>, ArrayGaps> >::Type & row1,
                Row<Align<String<Dna5>, ArrayGaps> >::Type & row2,
				String<Dna5> & genome,
                String<Dna5> & read, 
                String<Dna5> & comrev_read,
                uint64_t & cord,
                int block_size = _default_block_size_,
                int band = _default_block_size_ / 2,
                int local_flag = 1,
                Score<int> scheme = _default_scheme_
               )
{
    uint64_t cord_end = _DefaultCord.shift(cord, block_size, block_size);
    int score = align_cord (row1, row2, genome, read, comrev_read, cord, cord_end, band, band, local_flag, _default_scheme_);
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
    //std::cout << "ch" << row1 << "\n";
    //std::cout << "ch" << row2 << "\n";
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
int mergeAlignCheck_(Row<Align<String<Dna5>,ArrayGaps> >::Type & row11,
                 Row<Align<String<Dna5>,ArrayGaps> >::Type & row12,
                 Row<Align<String<Dna5>,ArrayGaps> >::Type & row21,
                 Row<Align<String<Dna5>,ArrayGaps> >::Type & row22,
                 uint64_t & cord1,  
                 uint64_t & cord2)
{
    int r_flag = 0;
    int64_t dx = get_cord_x(cord2) - get_cord_x(cord1);
    int64_t dy = get_cord_y(cord2) - get_cord_y(cord1);
    if (cord1 == emptyCord)
    {
        return 32;
    }
    int thd_cord_overlap = std::max(endPosition(row11), endPosition(row12));
    if (get_cord_strand(cord1 ^ cord2))
    {
        return  2;
    }
    //else if (!_DefaultCord.isCordsOverlap(cord1, cord2, thd_cord_overlap))  //sv: gap or reverse
    else if ((dx > thd_cord_overlap && dy > thd_cord_overlap)||
            (dx < 0 || dy < 0))
    {
        //print_cord(cord1, "man111");
        //print_cord(cord2, "man111");
        //dout << "man111" << thd_cord_overlap << get_cord_x(cord1) + endPosition(row11) << get_cord_x(cord2) + endPosition(row21) << "\n";
        //return 64 | 1;
    }
     
    int64_t delta1 = get_cord_x(cord2) - get_cord_x(cord1);
    int64_t delta2 = get_cord_y(cord2) - get_cord_y(cord1);
    if (int64_t(endPosition(row11) - beginPosition(row21)) < delta1 &&
        int64_t(endPosition(row12) - beginPosition(row22)) < delta2)
    {
        //print_cord(cord1, "manh");
        //print_cord(cord2, "manh");
//        return 1;
        //print_cord(cord1, "manh");
        //print_cord(cord2, "manh");
        //std::cout << "manh " << row11 << "\n";
        //std::cout << "manh " << row12 << "\n";
        //std::cout << "manh " << row21 << "\n";
        //std::cout << "manh " << row22 << "\n";
        //dout << "manh " << delta1 << endPosition(row11) << beginPosition(row21) << 
        //delta2 << endPosition(row12) << beginPosition(row22)  << "\n";
        //std::cout << "manh " << endPosition(row12) << " " << beginPosition(row22) + delta2 << "\n";
        return 1|4096;
    }
    if (int64_t(beginPosition(row11) - beginPosition(row21)) > delta1)
    {
        //return 1|128;
    }
    if (int64_t(beginPosition(row12) - beginPosition(row22)) > delta2)
    {
        //return 1|256;
    }
    if (int64_t(endPosition(row11) - endPosition(row21)) > delta1)
    {
        //return 1|1024;
    }
    if (int64_t(endPosition(row12) - endPosition(row22)) > delta2)
    {
        //return 1|2048;
    }
    return 0;
}

int mergeAlignCache1_(Row<Align<String<Dna5>,ArrayGaps> >::Type & row11,
                 Row<Align<String<Dna5>,ArrayGaps> >::Type & row12,
                 Row<Align<String<Dna5>,ArrayGaps> >::Type & row21,
                 Row<Align<String<Dna5>,ArrayGaps> >::Type & row22,
                 String<int64_t> & align1,
                 String<int64_t> & align2,
                 uint64_t & cord1,  
                 uint64_t & cord2)
{
    int bit = 20, bit2 = 40;
    int64_t mask = (1ULL << bit) - 1;
    int64_t delta1 = get_cord_x(cord2) - get_cord_x(cord1);
    int64_t delta2 = get_cord_y(cord2) - get_cord_y(cord1);
    TRowIterator it1, it2;
    //cache source coordinates of align row1
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
    
    //cache source coordinates of align row2
    it1 = begin(row21) + intersect_view_Begin;
    it2 = begin(row22) + intersect_view_Begin;
    for (int64_t i = intersect_view_Begin; i < intersect_view_End; i++)
    {
        if (*it1 == *it2)
        {
            appendValue (align2, ((i + clippedBeginPosition(row21))<< bit2) + 
                                 ((src1 + delta1) << bit) + (src2 + delta2));
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
        //dout << "mc12" << intersect_view_Begin << intersect_view_End << "\n";
        //std::cout << "mc1" << row11 << "\n";
        //std::cout << "mc1" << row12 << "\n";
        //std::cout << "mc1" << row21 << "\n";
        //std::cout << "mc1" << row22 << "\n";
        return 8 | 1;
    }
    return 0;
}
int mergeAlign1_(Row<Align<String<Dna5>,ArrayGaps> >::Type & row11,
                 Row<Align<String<Dna5>,ArrayGaps> >::Type & row12,
                 Row<Align<String<Dna5>,ArrayGaps> >::Type & row21,
                 Row<Align<String<Dna5>,ArrayGaps> >::Type & row22,
                 String<int64_t> & align1,
                 String<int64_t> & align2)
{
    int bit = 20, bit2 = 40;
    int flag = 0, start_j = 0;
    int64_t mask = (1ULL << bit) - 1;
    int64_t x1 = (align1[0] >> bit) & mask;
    int64_t y1 = align1[0] & mask;flag = 0;
    for (int i = 0; i < length(align1) - 1; i++)    
    {
        int64_t x1_next = (align1[i + 1] >> bit) & mask;
        int64_t y1_next = align1[i + 1] & mask;
        flag = 0;
        for (int j = start_j; j < length(align2); j++)  
        {
            int64_t x2 = align2[j] >> bit & mask;
            int64_t y2 = align2[j] & mask;
            int64_t dx = x2 - x1;
            int64_t dy = y2 - y1;
            if (!flag)
            {
                if (std::abs(x1_next - x2) < 2) 
                {
                    start_j = j;
                    flag = 1;
                }
            }
            //if (std::abs(x1 - x2) < thd_merge_x && std::abs(y1 - y2) < thd_merge_y)
            int n = 1;
            //if (dx + dy <= 1 && dx >=0 && dy >= 0)
            if (((dx == 0 && dy <= n) || (dx <= 6 && dy == 0)) && dx >= 0 && dy >=0)
            {
                int clip1 = (align1[i] >> bit2 & mask);
                int clip2 = (align2[j] >> bit2 & mask);
                //if (dy == 1)
                if (dx == 0)
                {
                    for (int i = 0; i < dy; i++)
                    {
                        insertGap(row11, clip1 - clippedBeginPosition(row11) - 1);
                        clip1++;
                    }
                }
                //else if (dx == 1)
                else if (dy == 0)
                {
                    for (int i = 0; i < dx; i++)
                    {
                        insertGap(row12, clip1 - clippedBeginPosition(row12) - 1);
                        clip1++;
                    }
                }
                setClippedBeginPosition(row21, clip2);
                setClippedBeginPosition(row22, clip2);
                setClippedEndPosition(row11, clip1);
                setClippedEndPosition(row12, clip1);
                //print_cord(cord1, "man21");
                //print_cord(cord2, "man22");
                return 0;
            }
            else if (x2 - x1 > 1 || y2 - y1 > 1)
            {
                break;
            }
        }
        x1 = x1_next;
        y1 = y1_next;
    }
    return 16 | 1;
}
/* 
 * Cache the coordinates(positions) to @align1 and @align2
   @align1:|empty[4]|c1[20]|c2[20]|c3[20]
 * c1[i] \in [0, length(row._array)) is unclipped view pos (not view or source pos). 
   It points to the row._array's c1[i]th bucket directly   
   use c1[i] - ClippedBeginPosition to get the corresponding view pos.
   c2[i] source pos of refs
   c3[i] source pos of reads.
 * @row1 and @row2 are supposed to be pair of rows of alignement, same size of view pos  
 * @delta1,@delta2 transform local src_x and src_y to global coordinates 
 */
int alignCachePos_(TRow5A & row1, TRow5A & row2, AlignCache & align,
                uint64_t src_x_l, uint64_t src_x_u, 
                uint64_t src_y_l, uint64_t src_y_u,
                int64_t cord_0)
{
    if (clippedBeginPosition(row1) != clippedBeginPosition(row2) ||
        clippedEndPosition(row1) != clippedEndPosition(row2))
    {
        return 1;
    }
    if (src_x_l > src_x_u && src_y_l > src_y_u)
    {
        return 1 | 2;
    }
    if (src_x_l < beginPosition(row1) || src_x_l > endPosition(row1) ||
        src_y_l < beginPosition(row2) || src_y_l > endPosition(row2))
    {
        return 1 | 4;
    }
    int bit = 20, bit2 = 40;
    RowPosViewer rpv1(row1), rpv2(row2);
    rpv1.findSrc(std::max(src_x_l, uint64_t(beginPosition(row1))));
    rpv2.findSrc(std::max(src_y_l, uint64_t(beginPosition(row2))));
    //dout << "capos2" << "view" << rpv1.getSrc()  << rpv1.getView() << rpv1.getUCView() << rpv1.getView() << "\n";
    if (rpv1.getView() < rpv2.getView())
    {
        rpv2.findView(rpv1.getView());
    }
    else if (rpv1.getView() > rpv2.getView())
    {
        rpv1.findView(rpv2.getView());
    }
    //dout << "capos3" << rpv2.getSrc() << rpv1.getView() << toSourcePosition(row2, rpv2.getView()) << rpv2.getView() << toSourcePosition(row1, rpv1.getView()) << "\n";
    while(rpv1.getSrc() < std::min(src_x_u, uint64_t(endPosition(row1))) || 
          rpv2.getSrc() < std::min(src_y_u, uint64_t(endPosition(row2))))
    {
        //rpv1.uc_pos1 == rpv2.up_view && uc_pos rpv2.view = rpv2.view
        align.appendValue(rpv1.getUCView(), rpv1.getSrc(), rpv2.getSrc(), cord_0);
        //appendValue (align, (rpv1.getUCView()<< bit2)
        //    + ((rpv1.getSrc() + delta1) << bit)
        //    + (rpv2.getSrc() + delta2));
        //dout << "capos1" << rpv1.getUCView() << rpv1.getView() << toViewPosition(*rpv1.rp, rpv1.getSrc())  << "\n";

        rpv1.nextView();
        rpv2.nextView();
    }
    if (align.empty())
    {
        return 8 | 1;
    }
    return 0;
}
int mergeAlignCache_(Row<Align<String<Dna5>,ArrayGaps> >::Type & row11,
                 Row<Align<String<Dna5>,ArrayGaps> >::Type & row12,
                 Row<Align<String<Dna5>,ArrayGaps> >::Type & row21,
                 Row<Align<String<Dna5>,ArrayGaps> >::Type & row22,
                 AlignCache & align1,
                 AlignCache & align2,
                 uint64_t & cord1,  
                 uint64_t & cord2)
{
    int bit = 20, bit2 = 40;
    int64_t mask = (1ULL << bit) - 1;
    uint64_t x1 = get_cord_x(cord1);
    uint64_t y1 = get_cord_y(cord1);
    uint64_t x2 = get_cord_x(cord2);
    uint64_t y2 = get_cord_y(cord2);
    uint64_t r1xl = x1 + beginPosition(row11);
    uint64_t r1xu = x1 + endPosition(row11);
    uint64_t r1yl = y1 + beginPosition(row12);
    uint64_t r1yu = y1 + endPosition(row12);
    uint64_t r2xl = x2 + beginPosition(row21);
    uint64_t r2xu = x2 + endPosition(row21);
    uint64_t r2yl = y2 + beginPosition(row22);
    uint64_t r2yu = y2 + endPosition(row22);
    uint64_t intersect_cord__xl = std::max(r1xl, r2xl); //row*1 lower bound (x)
    uint64_t intersect_cord__xu = std::min(r1xu, r2xu);
    uint64_t intersect_cord__yl = std::max(r1yl, r2yl); //row*2 lower bound (y)
    uint64_t intersect_cord__yu = std::min(r1yu, r2yu);


    if (intersect_cord__xl > intersect_cord__xu &&
        intersect_cord__yl > intersect_cord__yu)
    {
        return 1 ;
    }
    int f_a = 0;

    uint64_t intersect_src_1xl = std::min(intersect_cord__xl, r1xu) - x1;
    uint64_t intersect_src_1xu = std::max(intersect_cord__xu, r1xl) - x1;
    uint64_t intersect_src_1yl = std::min(intersect_cord__yl, r1yu) - y1;
    uint64_t intersect_src_1yu = std::max(intersect_cord__yu, r1yl) - y1;
    f_a |= alignCachePos_(row11, row12, align1, intersect_src_1xl, intersect_src_1xu, 
        intersect_src_1yl, intersect_src_1yu, cord1);

    uint64_t intersect_src_2xl = std::min(intersect_cord__xl, r2xu) - x2;
    uint64_t intersect_src_2xu = std::max(intersect_cord__xu, r2xl) - x2;
    uint64_t intersect_src_2yl = std::min(intersect_cord__yl, r2yu) - y2;
    uint64_t intersect_src_2yu = std::max(intersect_cord__yu, r2yl) - y2;
    f_a |= alignCachePos_(row21, row22, align2, intersect_src_2xl, intersect_src_2xu, 
        intersect_src_2yl, intersect_src_2yu, cord2);
    //dout << "macs" << intersect_src_2xl << intersect_src_2xu << intersect_src_2yl << intersect_src_2yu << beginPosition(row21) << "\n";
    //dout << "macs2" << f_a << "\n";
    return f_a; 
}
/**
 * Merge any size of two blocks which starts from cord1 and cord2.
 * @row_ij are supposed to contain the alignment of the two blocks.
 * Rows of @row_ij are compared and merged starting from clippedBeginPosition() 
   to clippedEndPosition() of each
 */
int merge_align_(Row<Align<String<Dna5>,ArrayGaps> >::Type & row11,
				 Row<Align<String<Dna5>,ArrayGaps> >::Type & row12,
				 Row<Align<String<Dna5>,ArrayGaps> >::Type & row21,
				 Row<Align<String<Dna5>,ArrayGaps> >::Type & row22,
                 String<Dna5> & ref,
                 String<Dna5> & read,
                 String<Dna5> & comrev_read,
				 uint64_t & cord1,  
				 uint64_t & cord2)
{
    //<<debug
    Row<Align<String<Dna5>, ArrayGaps> >::Type r1 = row11; 
    Row<Align<String<Dna5>, ArrayGaps> >::Type r2 = row12; 
    Row<Align<String<Dna5>, ArrayGaps> >::Type r3 = row21; 
    Row<Align<String<Dna5>, ArrayGaps> >::Type r4 = row22; 
    uint64_t c1 = cord1, c2 = cord2;
    //dout << "maa1" << beginPosition(r3) << beginPosition(r4) << "\n";
    //if (get_cord_y(cord2) == 2957 || get_cord_y(cord2) == 3255)
    //{
     //   print_cord(get_cord_y(cord2), "ins");
     //   std::cout << "ins11" << r1 << "\n";
     //   std::cout << "ins12" << r2 << "\n";
      //  std::cout << "ins21" << r3 << "\n";
       // std::cout << "ins22" << r4 << "\n";
    //}
    //print_cord(cord1, "man11");
    //print_cord(cord2, "man12");
    //>>debug
    int r_flag = 0;
    

    int f_c = mergeAlignCheck_(row11, row12, row21, row22, cord1, cord2);
    if (f_c)
    {
        //dout << "manx2<<" << f_c << get_cord_y(cord1) << get_cord_y(cord2) << "\n";
        //print_cord(cord1, "manx2");
        //print_cord(cord2, "manx2");

        return f_c;
    }

    //String<int64_t> align1, align2;
    AlignCache align1, align2; 
    f_c = mergeAlignCache_(row11, row12, row21, row22, align1, align2, cord1, cord2);
    if (f_c)
    {
        //dout << "manx3" << f_c << get_cord_y(cord1) << get_cord_y(cord2) << "\n";
        //print_cord(cord1, "manx3");
        //print_cord(cord2, "manx3");
        return f_c;
    }

    //<<debug
    //print_cord(cord1, "man60");
    //print_cord(cord2, "man61");
    //String<int64_t> a1, a2;
    //int f_ = mergeAlignCache_(r1, r2, r3, r4, a1, a2, cord1, cord2);
    //if (!f_)
    //{
    //    mergeAlign2_(r1, r2, r3, r4, a1, a2, ref, read, comrev_read, c1, c2);
    //}
    //>>debug
    //f_c = mergeAlign1_(row11, row12, row21, row22, align1, align2);
    f_c = mergeAlign2_(row11, row12, row21, row22, align1, align2, ref, read, comrev_read, cord1, cord2);
    //print_cord(cord1, "man11");
    //print_cord(cord2, "man12");
    
    //std::cout << "man13" << row11 << "\n";
    //std::cout << "man13" << row12 << "\n";
    //std::cout << "man13" << row21 << "\n";
    //std::cout << "man13" << row22 << "\n";
    //<<debug
        //printAlign_(row11, row12);
    //if (get_cord_y(c1) == 2957 || get_cord_y(c2) == 3255)
    //{

        //Align<String<Dna5>,ArrayGaps> align1; 
     //   std::cout << "m11" << r1 << "\n";
     //   std::cout << "m11" << r2 << "\n";
      //  std::cout << "m11" << r3 << "\n";
      //  std::cout << "m11" << r4 << "\n";

    //}
        //>>debug
    //dout << "manx1" << f_c << clippedBeginPosition(row11) - clippedBeginPosition(row12) << "\n";
    return f_c;
}

/**
 * if the new gap:=(@cord_str, @cord_end) can be merged with 
 * the last gap existing in the @gaps return 0;
 * else return1;
 * gap will joint with its adjacent if @f_merge is ture and interval between the two gaps < thd_merge_gap;
 */
int insertGaps(GapRecords & gaps,
               uint64_t cord_str, //start cord 
               uint64_t cord_end, //end cord 
               int bam_segs_id,
               int thd_merge_gap,
               int f_merge)
{
    int flag = 1;
    uint64_t cord1 = cord_str;
    uint64_t cord2 = cord_end;
    if (empty(gaps.c_pairs))
    {
        appendValue(gaps.c_pairs, std::pair<uint64_t, uint64_t>(cord1, cord2));
        appendValue(gaps.bam_segs_id, bam_segs_id);
    }
    else 
    {
        if (_DefaultCord.isCordsOverlap(back(gaps.c_pairs).second, cord1, thd_merge_gap) && f_merge)
        {
            int i = length(gaps.c_pairs) - 1;
            gaps.get_c_pair(i).second = cord2;
            flag = 0;
        }
        else
        {
            appendValue(gaps.c_pairs, std::pair<uint64_t, uint64_t>(cord1, cord2));
            appendValue(gaps.bam_segs_id, bam_segs_id);
        }
    }
    return flag;
}
int insertGaps(GapRecords & gaps,
               uint64_t cord_str, //start cord 
               uint64_t cord_end, //end cord 
               Row<Align<String<Dna5>, ArrayGaps> >::Type & row11,
               Row<Align<String<Dna5>, ArrayGaps> >::Type & row12,
               Row<Align<String<Dna5>, ArrayGaps> >::Type & row21,
               Row<Align<String<Dna5>, ArrayGaps> >::Type & row22,
               int bam_segs_id,
               int thd_merge_gap,
               int f_merge
                 )
{
    int flag = 1;
    uint64_t cord1 = cord_str;
    uint64_t cord2 = cord_end;
    if (empty(gaps.c_pairs))
    {
        appendValue(gaps.c_pairs, std::pair<uint64_t, uint64_t>(cord1, cord2));
        appendValue(gaps.r_pairs, std::pair<TRow, TRow>(row11, row12));
        appendValue(gaps.r_pairs, std::pair<TRow, TRow>(row21, row22));
        appendValue(gaps.bam_segs_id, bam_segs_id);
    }
    else 
    {
        if (_DefaultCord.isCordsOverlap(back(gaps.c_pairs).second, cord1, thd_merge_gap) && f_merge)
        {
            int i = length(gaps.c_pairs) - 1;
            gaps.get_c_pair(i).second = cord2;
            gaps.get_r2_pair(i).first = row21;
            gaps.get_r2_pair(i).second = row22;
            flag = 0;
        }
        else
        {
            appendValue(gaps.c_pairs, std::pair<uint64_t, uint64_t>(cord1, cord2));
            appendValue(gaps.r_pairs, std::pair<TRow, TRow>(row11, row12));
            appendValue(gaps.r_pairs, std::pair<TRow, TRow>(row21, row22));
            appendValue(gaps.bam_segs_id, bam_segs_id);
        }
    }
    return flag;
}

/*debug utility function*/
void printCigarSrcLen(String<BamAlignmentRecordLink> & records, CharString header = "pcsl ")
{
    for (int i = 0; i < length(records); i++)
    {
        std::pair <int, int > lens = countCigar (records[i].cigar);
        std::cout << header << " " << i << " " << length(records) << " " << lens.first << " " << lens.second << "\n";
    }
    std::cout << header << "\n";
}

/**
 *  Clip @row1 and @row2 into segments represented by pair of start and end
 *  of the segments. Each segment is supposed to be well aligned (of high align score).
 *  The size of the block to be clipped are supposed to > window_size,
 *  since the function clips the middle part of the block excluding 
 *  the head and tail part. Too small block can't be well clipped. 
 *  the result @clip_records @clip_records_src1 @clip_records_src2 are global 
 *  coordinates rather than relative coordinates in the rows
 */
int clip_rows_segs_(Row<Align<String<Dna5>,ArrayGaps> >::Type & row1,
                    Row<Align<String<Dna5>,ArrayGaps> >::Type & row2,
                    String<std::pair<int, int> > & clip_records,     //view coordinate
                    String<std::pair<int, int> > & clip_records_src1,//source1 coordinate
                    String<std::pair<int, int> > & clip_records_src2,//source2 coordinate
                    uint64_t gap_start_cord,
                    int thd_clip_score,
                    int thd_reject_score,
                    int thd_accept_score,
                    int thd_min_interval,
                    float thd_accept_density,
                    int thd_clip_mini_region = 20,
                    int window = 30, //supposed to < 50
                    int delta = 3)
{
    int x = 0;
    int thd_1 = 1;
    int flag = 0;
    clear(clip_records);
    clear(clip_records_src1);
    clear(clip_records_src2);
    int64_t src1 = beginPosition(row1);
    int64_t src2 = beginPosition(row2);
    TRowIterator it1 = begin(row1);
    TRowIterator it2 = begin(row2);
    TRowIterator it1_2 = it1, it2_2 = it2; 
    String<int> buffer;      //score at the intervals by delta
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
        if (k > 0 && !isGap(it1 - 1))
        {
            src1++;
        }
        if (k > 0 && !isGap(it2 - 1))
        {
            src2++;
        }
        if (count++ == delta)
        {
            buffer[buf_len] = x; //+ thd_1 * (std::abs(src1 - src1_2 - src2 + src2_2));
            buffer_src1[buf_len] = src1;
            buffer_src2[buf_len] = src2;
            //if (src1 != toSourcePosition(row1, buf_len * delta) ||
            //    src2 != toSourcePosition(row2, buf_len * delta))
            buf_len++;
            count = 1;
        }
        ++it1;
        ++it2;
    }
//TODO::change the score so it will not clip ins or dels.
    resize (buffer, buf_len);
    int prebp = 0;  //pre-breakpoint counter of buffer
    int bp = prebp; //bp * delta = view coordinate
    int d_w = window / delta;
    int d_m = thd_clip_mini_region / delta ;
//TODO!!::clip the fist and last segment
    int last_ = length(buffer) - std::max(d_w, d_m);
    int ct_clips = 0;  //count number of clip records
    int last_clip = 0;
    std::pair<int, int> clip_pair;
    std::pair<int, int> clip_pair_src1;
    std::pair<int, int> clip_pair_src2;
    for (int i = d_w + 1; i < last_ ; i++)
    {
        int d_score_left = buffer[i] - buffer[i - d_w];
        int d_score_right = buffer[i + d_w] - buffer[i];
        int d_score = buffer[i + d_w] - (buffer[i] << 1) + buffer[i - d_w];
        if ((std::abs(d_score) > thd_clip_score && 
             std::min(d_score_left, d_score_right) < thd_reject_score)|| 
             i == last_ - 1) 
        {
            bp = i;
            int max_d_score = std::abs(d_score);
            int d_j = std::min(i + d_m, (int)length(buffer) - d_w);
            for (int j = i + 1; j < d_j; j++)
            {
                d_score = buffer[j + d_w] - (buffer[j] * 2) + buffer[j - d_w]; 
                if (std::abs(d_score) > max_d_score)
                {
                    bp = j;
                    max_d_score = std::abs(d_score);
                } 
            } //retrieve the max clip_score within the thd_clip_mini_region 
            float dens_ = (float)(buffer[bp] - buffer[prebp]) / 
                                 (delta * (bp - prebp + 1));
//TODO::score to accecpt the segs, need parm tunning 
            if (dens_ > thd_accept_density)
            {
                clip_pair.first = prebp * delta + clipped_begin;
                clip_pair.second = bp * delta + clipped_begin;
                if (empty(clip_records) || 
                    clip_pair.first - last_clip > thd_min_interval)
                {
                    clip_pair_src1.first = get_cord_x(gap_start_cord) + buffer_src1[prebp];
                    clip_pair_src1.second = get_cord_x(gap_start_cord) + buffer_src1[bp];
                    clip_pair_src2.first = get_cord_y(gap_start_cord) + buffer_src2[prebp];
                    clip_pair_src2.second = get_cord_y(gap_start_cord) + buffer_src2[bp];
                    appendValue(clip_records, clip_pair);
                    appendValue(clip_records_src1, clip_pair_src1);
                    appendValue(clip_records_src2, clip_pair_src2);
                    ++ct_clips;
                }
                else
                {
                    back(clip_records).second = clip_pair.second;
                    back(clip_records_src1).second = clip_pair_src1.second;
                    back(clip_records_src2).second = clip_pair_src2.second;
                }
                last_clip = clip_pair.second; 
            }
            prebp = bp;
            i = d_j;
        }
    }

    return 0;
}
int clip_rows_segs (Row<Align<String<Dna5>,ArrayGaps> >::Type & row1,
               Row<Align<String<Dna5>,ArrayGaps> >::Type & row2,
               String<std::pair<int, int> > & clip_records,     //view coordinate
               String<std::pair<int, int> > & clip_records_src1,//source1 coordinate
               String<std::pair<int, int> > & clip_records_src2,//source2 coordinate
               uint64_t gap_start_cord,
               AlignGapParms align_gap_parms,
               int thd_clip_mini_region = 20,
               int window = 30, //supposed to < 50
               int delta = 3
              )
{
    return clip_rows_segs_(row1, row2,
                    clip_records,     //view coordinate
                    clip_records_src1,//source1 coordinate
                    clip_records_src2,//source2 coordinate
                    gap_start_cord,
                    align_gap_parms.thd_clip_score,
                    align_gap_parms.thd_reject_score,
                    align_gap_parms.thd_accept_score,
                    align_gap_parms.thd_min_interval,
                    align_gap_parms.thd_accept_density,
                    thd_clip_mini_region,
                    window, //supposed to < 50
                    delta
                   );
}
/*
 * Clip head or tail by calling the clip_rows_segs
 * @direction -1 head, 1 tail
 */
int clip_segs(Row<Align<String<Dna5>,ArrayGaps> >::Type & row1,
               Row<Align<String<Dna5>,ArrayGaps> >::Type & row2,
               uint64_t gap_start_cord,
               AlignGapParms & align_gap_parms,
               int direction,   
               int thd_clip_mini_region = 20,
               int window = 30, //supposed to < 50
               int delta = 3
              )
{
    String<std::pair<int, int> > clip_records;
    String<std::pair<int, int> > clip_records_src1;
    String<std::pair<int, int> > clip_records_src2;
    clip_rows_segs_(row1, 
                   row2, 
                   clip_records, 
                   clip_records_src1, 
                   clip_records_src2, 
                   gap_start_cord, 
                   align_gap_parms.thd_clip_score,
                   align_gap_parms.thd_reject_score, 
                   align_gap_parms.thd_accept_score, 
                   align_gap_parms.thd_min_interval,
                   align_gap_parms.thd_accept_density
                );
    if (empty(clip_records))
    {
        return 1;
    }
    else 
    {
        if (direction <= 0)
        {
            setClippedBeginPositions(row1, row2, clip_records[0].first);
        }
        else if (direction >= 0)
        {
            setClippedEndPositions(row1, row2, back(clip_records).second);
        }
    }
    return 0;
}
/*-------  To clip head and tail of records in String<BamRecordLink>  -------*/
struct AlignClipScores
{
    int precision;
    int thd_window_size;
    int thd_edge_window_size;
    int thd_dens_match_lower;
    int thd_dens_match_upper;
    int thd_ddens_lower;
    String<int> scores; 
    String<int> views;
    String<int> src1;
    String<int> src2;
    String<int> matches;

    AlignClipScores();
    int & operator [](int);
    int cigar2Score(CigarElement<> & cigar);
    int appendNew(CigarElement<> & cigar);
    int addNew(CigarElement<> & cigar, int i);
    int findSrc2(int src2_value);
    void resize(int len);
};
AlignClipScores::AlignClipScores():
    precision(1000),
    thd_window_size (50),
    thd_edge_window_size(10),
    thd_dens_match_lower(0.65 * precision),
    thd_dens_match_upper(0.75 * precision),
    //the minium ddens that will be accepted
    //suppose lower side >40% gaps, upper sider < 20% gaps;
    //-10 per gaps and 1 per match in average, 
    //thd_ddens_lower= -(0.4*(-10)+0.6) + (0.2*(-10)+0.8))  = 2.2
    thd_ddens_lower(2.2 * precision) //relate to function cigar2Score
{}
int & AlignClipScores::operator [](int i)
{
    return scores[i];
}
/*
 * 'M' of cigar not allowed, '=' and 'X' required.
 */
int AlignClipScores::cigar2Score(CigarElement<> & cigar)
{
    int score = 0;
    if (cigar.operation == 'X')
    {
        score = -15 - cigar.count;
    }
    else if (cigar.operation == 'I' || cigar.operation == 'D')
    {
        if (cigar.count < 20)
        {
            score = -10 - cigar.count / 5;
        }
        else
        {
            score = -10 - cigar.count / 10;
        }
    }
    else if (cigar.operation == '=')
    {
        score = cigar.count;
    }
    else 
    {
        score = 0;
    }
    return score * precision; 
}
int AlignClipScores::appendNew(CigarElement<> & cigar)
{
    //dout << "acs" << "\n";
    int new_score, new_match, new_view;
    if (empty(scores))
    {
        new_score = new_view = new_match = 0;
    }
    else
    {
        new_score = back(scores);
        new_view = back(views);
        new_match = back(matches);
    }
    new_score += cigar2Score(cigar);
    new_view += cigar.count;
    if (cigar.operation == '=')
    {
        new_match += cigar.count * precision;
    }
    appendValue(scores, new_score);
    appendValue(views, new_view);
    appendValue(matches, new_match);
    //dout << "acs1" << back(views) << "\n";
    //dout << "appendnew" << back(scores) << back(views) << back(matches) << length(scores) << "\n";
    return 0;
}
int AlignClipScores::addNew(CigarElement<> & cigar, int i)
{
    //dout << "acs" << "\n";
    int new_score, new_match, new_view;
    if (empty(scores) || i == 0)
    {
        new_score = new_view = new_match = 0;
    }
    else
    {
        new_score = scores[i - 1];
        new_view = views[i - 1];
        new_match = matches[i - 1];
    }
    new_score += cigar2Score(cigar);
    new_view += cigar.count;
    if (cigar.operation == '=')
    {
        new_match += cigar.count * precision;
    }
    scores[i] = new_score;
    views[i] = new_view;
    matches[i] = new_match;
    //dout << "addnew" << scores[i] << views[i] << matches[i] << i + 1<< "\n";
    return 0;
}
int AlignClipScores::findSrc2(int src2_value)
{
    for (int i = 0; i < length(src2); i++)
    {
        if (src2[i] >= src2_value)
        {
            return i;
        }
    }
    return length(src2) - 1;
}
void AlignClipScores::resize(int len)
{
    seqan::resize (scores, len);
    seqan::resize (src2, len);
    seqan::resize (views, len);
    seqan::resize (matches, len);
}
/*
 * Struct called by internal functions only to record ranges of read for clipping
 * [@head_str, @head_end)... of read  are supposed to be the region to be clipped
 */
struct _HeadTailRange
{
    unsigned it; // itth head of BamRecordsLink
    int head_str;
    int head_end;
    int tail_str;
    int tail_end;

    void setStrEnd(int y_str1, int y_end1, int y_str2, int y_end2);
    void revertStrEnd(int read_len); //revert to different strand
    _HeadTailRange();
    _HeadTailRange(unsigned i, int y_str, int y_end);
};
void _HeadTailRange::setStrEnd(int y_str1, int y_end1, int y_str2, int y_end2)
{
    head_str = y_str1;
    head_end = y_end1;
    tail_str = y_str2;
    tail_end = y_end2;
}
void _HeadTailRange::revertStrEnd(int read_len)
{
    int tmp_str = read_len - 1 - tail_end;
    int tmp_end = read_len - 1 - tail_str;
    tail_str = read_len - 1 - head_end;
    tail_end = read_len - 1 - head_str;
    head_str = tmp_str; 
    head_end = tmp_end;
}
_HeadTailRange::_HeadTailRange(){}
_HeadTailRange::_HeadTailRange(unsigned i, int y_str, int y_end)
{
    setStrEnd(y_str, y_end, y_str, y_end);
    it = i;
}
int _bamRecordLlink2Score(String<BamAlignmentRecordLink> & records,
                          AlignClipScores & scores,
                          _HeadTailRange & range,
                          int it)
{
    //dout << "brl2s<<<" << "\n";
    int new_score = 0;
    //int new_src1 = records[it].beginPos;
    int seg_str = 0;
    int seg_len = 0;
    int new_view = 0;
    int it_origin = it;
    int len = 0;
    
    while(true)
    {
        //printCigars(records[it].cigar, "brcigar1");
        len += length(records[it].cigar);
        if (records[it].isEnd())
        {
            break;
        }
        else
        {
            it = records[it].next();
        }
    }
    scores.resize(len);
    
    it = it_origin;
    int ii = 0;
    int f_first_S = 1;
    while(true)
    {
        for (int i = 0; i < length(records[it].cigar); i++)
        {
            scores.addNew(records[it].cigar[i], ii);
            char o = records[it].cigar[i].operation;
            if ((o == 'S' || o == 'H') && f_first_S)
            {
                seg_str += records[it].cigar[i].count;
                f_first_S = 0;
            }
            else if (o == 'X' || o == 'I' || o == '=')
            {
                seg_len += records[it].cigar[i].count;
            }
            scores.src2[ii] = seg_str + seg_len;
            ++ii;
        }
        if (records[it].isEnd())
        {
            break;
        }
        else
        {
            it = records[it].next();
        }
    }
    //range = _HeadTailRange(seg_str, seg_str + seg_len + 1);
    range.setStrEnd(seg_str, seg_str + seg_len + 1, seg_str, seg_str + seg_len + 1);
    return 0;
}
/*
 * @f_ht<0 clip head, >0 to clip tail
 */
int _clipAlignScore(AlignClipScores & scores,
                    int view_str, 
                    int view_end,
                    int f_ht,
                    AlignGapParms & align_gap_parms)
{
    //dout << "csc" << length(scores.scores) << "\n";
    if (length(scores.scores) < 3 || view_str >= view_end)
    {
        return f_ht < 0 ? view_str : view_end; //error and do nothing
    }
    int clip = 0;
    if (f_ht < 0)
    {
        f_ht = -1;
        clip = view_str;
    }
    else if (f_ht > 0)
    {
        f_ht = 1;
        clip = view_end - 1;
    }
    int max_d_dens = scores.thd_ddens_lower;
    int i0 = view_str, i1 = view_str, i2 = view_str;
    for (i1 = view_str; i1 < view_end; i1++)
    {
        //dout << "cbx2" << scores.views[i0] << scores.views[i1] << scores.views[i2] << scores.thd_window_size << "\n";
        int f1 = 0, f2 = 0;
        while(i0 < i1 && scores.views[i1] - scores.views[i0] >= scores.thd_edge_window_size)
        {
            if ((scores.views[i1] - scores.views[i0] >= scores.thd_window_size &&
                 scores.views[i1] - scores.views[i0 + 1] < scores.thd_window_size) ||
                (scores.views[i1] - scores.views[i0] <  scores.thd_window_size && i0 == view_str))
            {
                f1 = 1;
                break;
            }
            i0++;
        }
        while (i2 < view_end)
        {
            if (scores.views[i2] - scores.views[i1] >= scores.thd_window_size || 
               (scores.views[i2] - scores.views[i1] >  scores.thd_edge_window_size && 
                scores.views[i2] - scores.views[i1] <  scores.thd_window_size && i2 == view_end - 1))
            {
                f2 = 1;
                break;
            }
            i2++;
        }
        if (f1 && f2)
        {
            int dens_left = (scores.scores[i1] - scores.scores[i0]) / 
                            (scores.views[i1] - scores.views[i0]);
            int dens_rght = (scores.scores[i2] - scores.scores[i1]) / 
                            (scores.views[i2] - scores.views[i1]);
            int dense_match_left = (scores.matches[i1] - scores.matches[i0]) / (scores.views[i1] - scores.views[i0]);
            int dense_match_rght = (scores.matches[i2] - scores.matches[i1]) / (scores.views[i2] - scores.views[i1]);
            //dout << "cbx1" << scores.views[i0] << scores.views[i1] << scores.views[i2] << i1 << (dens_left - dens_rght) * f_ht << dens_left << dens_rght << scores.scores[i1] - scores.scores[i0] << scores.scores[i2] - scores.scores[i1] << src1[i1] << src2[i1] << "match" <<scores.matches[i1] - scores.matches[i0] << scores.matches[view_end - 1] - scores.matches[i1] << "left" << dense_match_left << dense_match_rght << scores.thd_dens_match_lower << scores.thd_dens_match_upper << i1 << clip << "\n";

            if (f_ht < 0)
            {
                int ddens = dens_rght - dens_left;
                if (ddens > max_d_dens &&
                    dense_match_left < scores.thd_dens_match_lower &&
                    dense_match_rght > scores.thd_dens_match_upper)
                {
                    max_d_dens = ddens;
                    clip = i1;
                }
            }
            else if (f_ht > 0)
            {
                int ddens = dens_left - dens_rght;
                if (ddens > max_d_dens && 
                    dense_match_rght < scores.thd_dens_match_upper &&
                    dense_match_left > scores.thd_dens_match_lower)
                {
                    max_d_dens = ddens;
                    clip = i1 + 1;
                }
            }
        }
    }
    //dout << "cbx3" << clip << src1[clip] << src2[clip] << view_end << "\n";
    return clip;
}
/*
 * Operator function to manipulate(erase) head of the @records[@it] 
   given the head position @cigar_clip
   region [head, @cigar_clip) is clipped
 */
int _clipBamRecordLinkCigarHead(String<BamAlignmentRecordLink> & records,
                                int cigar_clip,
                                int it)
{
    int l = 0;
    std::pair<int, int> seqs_len = std::pair<int, int>(0, 0);
    std::pair<int, int> tmp;
    int first_it = it;
    int read_str = 0;
    int cigar_erased = 0;
    char read_opt = 'S';
    if (!empty(records[it].cigar))
    {
        CigarElement<> & cigar0 = records[it].cigar[0];
        if (cigar0.operation == 'S' || cigar0.operation == 'H')
        {
            read_str = cigar0.count;
            read_opt = cigar0.operation;
        }
    }
    while (true)
    {
        l += length(records[it].cigar);
        if (l > cigar_clip)
        {
            tmp = cigars2SeqsLen(records[it].cigar, 0, cigar_clip - l + length(records[it].cigar));
            seqs_len.first += tmp.first;
            seqs_len.second += tmp.second;
            //std::cout << records[it].cigar[]
            //dout << "cbx6" << cigar_clip - l + length(records[it].cigar) << "\n";
            if (cigar_clip - l + length(records[it].cigar) <= 1)
            {
                if (records[it].cigar[0].operation == 'S' || records[it].cigar[0].operation == 'H')
                {
                    ////std::cout << "cbx5" << read_opt << " " << read_str + seqs_len.second << "\n";
                    records[it].cigar[0].count = read_str + seqs_len.second;
                }
                else if (read_str + seqs_len.second != 0)
                {
                    //std::cout << "cbx4" << read_opt << " " << read_str + seqs_len.second << "\n";
                    insertValue(records[it].cigar, 0, CigarElement<>(read_opt, read_str + seqs_len.second));
                    cigar_erased += -1;
                }
            }
            else if (cigar_clip - l + length(records[it].cigar) > 1)
            {
                cigar_erased += cigar_clip - l + length(records[it].cigar) - 1;
                erase(records[it].cigar, 1, cigar_clip - l + length(records[it].cigar));
                records[it].cigar[0] = CigarElement<>(read_opt, read_str + seqs_len.second);
            }

            //dout << "cliphead" << cigar_clip << l << length(records[it].cigar) << cigar_clip - l + length(records[it].cigar) << "\n";
            //insertValue(records[it].cigar, 0, CigarElement<>(read_opt, read_str + seqs_len.second));
            break;
        }
        else
        {
            tmp = cigars2SeqsLen(records[it].cigar, 0, length(records[it].cigar));
            seqs_len.first += tmp.first;
            seqs_len.second += tmp.second;
        }
        if (records[it].isEnd())
        {
            break;
        }
        else
        {
            cigar_erased += length(records[it].cigar);
            clear(records[it].cigar);
            it = records[it].next();
        }
    }
    records[first_it].beginPos += seqs_len.first;
    //dout << "chhh2" << seqs_len.first << records[first_it].beginPos << cigar_erased << "\n";
    return cigar_erased;
}
/*
 * Operator function to manipulate(erase) tail of the @records[@it] 
   given the tail position @cigar_clip
   region [@cigar_clip, end) is clipped
 */
int _clipBamRecordLinkCigarTail(String<BamAlignmentRecordLink> & records,
                                int cigar_clip,
                                int it)
{
    int f_clear = 0;
    int l = 0;
    //dout << "cct1" << cigar_clip << "\n";
    while (true)
    {
        //dout << "cct3" << l << cigar_clip << length(records[it].cigar)<< "\n";
        l += length(records[it].cigar);
        if (f_clear)
        {
            clear(records[it].cigar);
        }
        else if (l >= cigar_clip)
        {
            //dout << "cct2" << cigar_clip << cigar_clip - l + length(records[it].cigar) << length(records[it].cigar)<< "\n";
            erase(records[it].cigar, cigar_clip - l + length(records[it].cigar), 
                length(records[it].cigar));
            f_clear = 1;
        }
        if (records[it].isEnd())
        {
            break;
        }
        else
        {
            it = records[it].next();
        }
    }
    return 0;
}
/*
 * clip head in region [str, str+@src_d_bps) and tail [end - @src_d_bps, end).
 * @src_d_bps short for source shift d base pairs, it's the length of bps going to be 
   clipped
 */
int _clipBamRecordLinkHeadTail(String<BamAlignmentRecordLink> & records, 
                               AlignClipScores & scores,
                               _HeadTailRange & range,
                               int it,
                               AlignGapParms & align_gap_parms)
{
    if (length(scores.scores) != length(scores.views) || empty(scores.scores))
    {
        return -1; //error 
    }
    int cigar_erased = 0;
    int clip_head = 0;
    for (int i = 0; i < length(scores.scores); i++)
    {
        if (scores.views[i] - scores.views[0] >= align_gap_parms.thd_clip_view_len  || 
            i == length(scores.scores) - 1)
        {
            int head_end_i = std::min(scores.findSrc2(range.head_end), i);
            int clip = _clipAlignScore(scores, 0, head_end_i, -1, align_gap_parms);
            cigar_erased = _clipBamRecordLinkCigarHead(records, clip, it);
            clip_head = clip;
            break;
        }
    }
    //dout << "chh1" << length(scores.scores) << clip_head << cigar_erased << "\n";
    for (int i = length(scores.scores) - 1; i > clip_head; i--)
    {
        if (back(scores.views) - scores.views[i] >= align_gap_parms.thd_clip_view_len || 
            i == clip_head + 1)
        {
            int tail_str_i = std::max(scores.findSrc2(range.tail_str), i);
            int clip = _clipAlignScore(scores, tail_str_i, length(scores.views), 1, align_gap_parms) - cigar_erased;
            _clipBamRecordLinkCigarTail(records, clip, it);
            //dout << "cht1" << clip << "\n";
            break;                
        }
    }
    return 0;
}
/*
 * Simple structs for recording value of sweeping and coverage.
   Supposed to be called locally by the function _calClipBamRecordLinkRange
 */
struct _SweepNode
{
    int key;
    int value;
    int type;

    _SweepNode();
    _SweepNode(int, int, int);
};
_SweepNode::_SweepNode(){}
_SweepNode::_SweepNode(int k, int v, int t) :
    key(k),
    value(v),
    type(t)
{}
struct _RangeCoverage
{
    int r_str;
    int r_end; 
    int coverage;

    _RangeCoverage();
    _RangeCoverage(int, int, int);
};
_RangeCoverage::_RangeCoverage(){}
_RangeCoverage::_RangeCoverage(int r1, int r2, int c) :
    r_str(r1),
    r_end(r2),
    coverage(c)
{}
/*
 * Calculate the range of clipping headd and tail for each line in the @records
 * Principle to determine the range: the region whose y (read) is not overlapped 
   with any other region of y, namely the coverage <=1, is not allowed to be 
   clipped.
 * @clip_ranges[i].first is the end of clipping head,
   @clip_ranges[i].second is the start of clipping tail.
 * Method: @range is supposed to contain the start and end of each segment 
   constituting the one complete record ind @records, the line. @sweep_nodes is 
   then initiated with  value ofstart and end from the @range. The @sweep_nodes 
   is then sorted and swept. If the value is value of start then coverage+=1,
   else if the value is value of end the coverage-=1. We then get coverages of 
   all subsegments. Then we determin the bound of x range, within which the 
   coverage >1, for clipping.
 */
int _calClipBamRecordLinkRange(String<BamAlignmentRecordLink> & records,
                               String<_HeadTailRange> & ranges,
                               unsigned read_len)
{
    unsigned thd_min_coverage_len = 20; // segment > 20bps of will be recorded
    BamLinkStringOperator bs;
    if (empty(ranges) || bs.getHeadNum(records) == 0)
    {
        return 0;
    }
    String<_SweepNode> sweep_nodes;
    resize (sweep_nodes, bs.getHeadNum(records) * 2);
    
    for (int i = 0; i < length(ranges); i++) //initiate with str and end of line
    {
        if (bs.getLineStrand(records, ranges[i].it))
        {
            ranges[i].revertStrEnd(read_len);
        }
        sweep_nodes[i * 2] = _SweepNode(i, ranges[i].head_str, 0);
        sweep_nodes[i * 2 + 1] = _SweepNode(i, ranges[i].tail_end, 1);
    }
    std::sort (begin(sweep_nodes), end(sweep_nodes), [](_SweepNode & a, _SweepNode & b){
            return a.value < b.value;
        });
    String<_RangeCoverage> ranges_coverage;
    unsigned coverage = 0; //coverage of given range of reads
    appendValue(ranges_coverage, _RangeCoverage(0, 0, 0));
    for (int i = 0; i < length(sweep_nodes); i++)
    {
        back(ranges_coverage).r_end = sweep_nodes[i].value;
        back(ranges_coverage).coverage = coverage;
        appendValue (ranges_coverage, _RangeCoverage(sweep_nodes[i].value, 0, 0));
        coverage = sweep_nodes[i].type == 0 ? coverage + 1 : coverage - 1;
    }
    erase(ranges_coverage, 0);  //the first and last are meaningless, thus erased
    erase(ranges_coverage, length(ranges_coverage) - 1);
    unsigned it = 0;
    for (int i = 0; i < length(ranges_coverage); i++)
    {
        if (ranges_coverage[i].coverage < 2 &&
            ranges_coverage[i].r_end - ranges_coverage[i].r_str > thd_min_coverage_len)
        {
            ranges_coverage[it] = ranges_coverage[i];
            if (it > 0 && 
                ranges_coverage[it].r_str - ranges_coverage[it - 1].r_end < thd_min_coverage_len)
            {
                ranges_coverage[it - 1].r_end = ranges_coverage[it].r_end;
                it--;
            }
            it++;
        }
    }
    //dout << "ccb1A" << it << "\n";
    resize(ranges_coverage, it);
    for (int i = 0; i < length(ranges); i++)
    {
        for (int j = 0; j < length(ranges_coverage); j++)
        {
            unsigned lower_bound = std::max(ranges[i].head_str, ranges_coverage[j].r_str);
            unsigned upper_bound = std::min(ranges[i].tail_end, ranges_coverage[j].r_end);
            //dout << "ccbr6" << j << ranges[i].head_str << ranges_coverage[j].r_str << lower_bound << upper_bound << "\n";
            if (lower_bound < upper_bound)
            {
                ranges[i].head_end = lower_bound;
                break;
            }
        }
    }
    for (int i = 0; i < length(ranges); i++)
    {
        for (int j = length(ranges_coverage) - 1; j >= 0; j--)
        {
            unsigned lower_bound = std::max(ranges[i].head_str, ranges_coverage[j].r_str);
            unsigned upper_bound = std::min(ranges[i].tail_end, ranges_coverage[j].r_end);
            if (lower_bound < upper_bound)
            {
                ranges[i].tail_str = upper_bound;
                break; 
            }
        } 
    }
    for (int i = 0; i < length(ranges); i++)
    {
        if (bs.getLineStrand(records, ranges[i].it)) //strand of it th line 
        {
            ranges[i].revertStrEnd(read_len);
        }
        //dout << "ccbr1" << bs.getLineStrand(records, ranges[i].it) << ranges[i].it << ranges[i].head_str << ranges[i].head_end << ranges[i].tail_str << ranges[i].tail_end << "\n";
    }
    return 0;
}
int clipBamRecordLinksHeadTail(String<BamAlignmentRecordLink> & records,
                               unsigned read_len,
                               AlignGapParms & align_gap_parms)
{
    BamLinkStringOperator bs;
    String<AlignClipScores> scores_set;
    String<_HeadTailRange> ranges;
    //StringSet<String<int> > src1s; 
    //StringSet<String<int> > src2s; 

    bs.updateHeadsTable(records);
    resize(ranges, bs.getHeadNum(records));
    resize(scores_set, bs.getHeadNum(records));
    //resize(src1s, bs.getHeadNum(records));
    //resize(src2s, bs.getHeadNum(records));
    for (unsigned i = 0; i < bs.getHeadNum(records); i++)
    {
        ranges[i] = _HeadTailRange(i, 0, 0);
        _bamRecordLlink2Score(records, scores_set[i], ranges[i], bs.getHead(records, i));
    }
    _calClipBamRecordLinkRange(records, ranges, read_len);
    for (int i = 0; i < bs.getHeadNum(records); i++)
    {
        _clipBamRecordLinkHeadTail(records, scores_set[i], ranges[i], 
            bs.getHead(records, ranges[i].it), align_gap_parms);
    }
    return 0;
}
/*----------  End  ----------*/
/*
 *!Make sure head and tail in the gap specified by @gap are well aligned.
 * Head and tail will be joined to the origin directly
 * such that the gap can be merged into the main alignment;
 *@thd_alg_extnd: extension in inversion sub region to be aligned
 * such that it can be clipped precisely at two ends;
 */
int align_gap (GapRecordHolder & gap,
               String<BamAlignmentRecordLink> & bam_records,
               StringSet<String<Dna5> >& genomes,
               String<Dna5> & read, 
               String<Dna5> & comrev_read,
               Score<int> & score_scheme,
               AlignGapParms & align_gap_parms)
{
    typedef String<std::pair<int, int> > ClipRecords;
    int thd_alg_extnd = 20;
    uint64_t str_cord = gap.getCords().first;
    uint64_t end_cord = gap.getCords().second;
    //return 0;
    if (get_cord_strand(str_cord ^ end_cord))
    {
        return 1;
    }
    ClipRecords clips, clips_src1, clips_src2;
    ClipRecords seg_clips, seg_clips_src1, seg_clips_src2;
    uint64_t seg_str_cord;
    uint64_t seg_end_cord;
    int g_id = get_cord_id(str_cord);
    int band = std::max(get_cord_x(end_cord) - get_cord_x(str_cord),
                        get_cord_y(end_cord) - get_cord_y(str_cord)) / 2;
    int bam_id = gap.getBamSegIdHead();
    int bam_next_id = gap.getBamSegIdTail();
    //WARNING::modify band::too large band
    TRow row1, row2, row3, row4 ;
    //dout << "ag3" << get_cord_id(str_cord) << get_cord_id(end_cord) << get_cord_y(str_cord) << get_cord_y(end_cord) << get_cord_x(str_cord) << get_cord_x(end_cord) << band << "\n";
    //return 0;
    align_cord (row1, row2, genomes[g_id], read, comrev_read, str_cord, end_cord, band, band);
    //Head and tail are already merged, so view_str = 0.
    int const view_str = 0;
    int const view_end = clippedEndPosition(row1);
    clip_rows_segs (row1, row2, clips, clips_src1, clips_src2, str_cord, align_gap_parms);
    //std::cout << "ag11" << row1 << "\n";
    //std::cout << "ag12" << row2 << "\n";
    //dout << "ag2" << length(clips) << get_cord_y(str_cord) << get_cord_y(end_cord) << "\n";
    //std::cout << "ag2" << row1 << "\n";
    //std::cout << "ag2" << row2 << "\n";
    if (empty(clips))
    {
        return 1;
    }
    if (length (clips) == 1)
    {
        setClippedPositions(row1, row2, view_str, view_end);
        insertBamRecordCigar(bam_records[bam_id], row1, row2);
        addNextBamLink (bam_records, bam_id, bam_next_id);
    }
    else if (length(clips) > 1)
    {
        setClippedPositions(row1, row2, view_str, clips[0].second);
        insertBamRecordCigar(bam_records[bam_id], row1, row2);
        setClippedPositions(row1, row2, back(clips).first, view_end);
        int bam_start_x = get_cord_x (str_cord) + beginPosition(row1);
        int bam_start_y = get_cord_y (str_cord) + beginPosition(row2);
        insertBamRecord(bam_records[bam_next_id], row1, row2, g_id, bam_start_x, bam_start_y, 0);
        //int tmp = bam_next_id;
        //int tmp2 = tmp + 1;
        int pre_bam_id = bam_id;
        for (int i = 1; i < length(clips) - 1; i++)
        {
            int bam_start_x = clips_src1[i].first;
            int bam_start_y = clips_src2[i].first;
            int bam_strand = get_cord_strand(str_cord);
            setClippedPositions(row1, row2, clips[i].first, clips[i].second);
            insertNewBamRecord(bam_records, row1, row2, g_id, bam_start_x, bam_start_y, bam_strand, -1, 1, 2048); 
            //addNextBamLink(bam_records, pre_bam_id, length(bam_records) - 1);
            //dout << "ag14" << bam_st<<  "\n";
            //std::cerr << "it=====" << pre_bam_id << length(bam_records) - 1 << "\n";
            pre_bam_id = i;
        }
        //addNextBamLink(bam_records, length(bam_records) - 1, bam_next_id);
        //std::cerr << "it=====" << length(bam_records) - 1 << bam_next_id << "\n";
        //realign and clip interval between each clip[i].second and clip[i + 1].second
        for (int i = 0; i < length(clips) - 1; i++) 
        {
            seg_str_cord = new_xy_cord (str_cord, 
                                        clips_src1[i].second, 
                                        clips_src2[i].second);
            seg_end_cord = new_xy_cord (end_cord, 
                                        clips_src1[i + 1].first, 
                                        clips_src2[i + 1].first);
            cmpRevCord (seg_str_cord, seg_end_cord, 
                        seg_str_cord, seg_end_cord, length(read));
            seg_str_cord = shift_cord(seg_str_cord, -thd_alg_extnd, -thd_alg_extnd);
            seg_end_cord = shift_cord(seg_end_cord, thd_alg_extnd, thd_alg_extnd);
            int seg_band = std::max(get_cord_x(seg_end_cord - seg_str_cord),
                                    get_cord_y(seg_end_cord - seg_str_cord)) / 2;
            align_cord (row3, row4, genomes[g_id], 
                        read, comrev_read, seg_str_cord, seg_end_cord, seg_band, seg_band);
            clip_rows_segs (row3, row4, 
                            seg_clips, 
                            seg_clips_src1, 
                            seg_clips_src2,
                            seg_str_cord,align_gap_parms);
            //int tmp3 = tmp2;
            //>>debug
            
            if (empty(seg_clips))
            {

            //dout << "ag18" << length(seg_clips) << "\n";
            int bam_start_x = clips_src1[i].second;
            int bam_start_y = clips_src2[i].second;
            int bam_strand = get_cord_strand(str_cord);
                setClippedPositions(row1, row2, clips[i].second, clips[i + 1].first);
            insertNewBamRecord(bam_records, row1, row2, g_id, bam_start_x, bam_start_y, bam_strand, -1, 1, 2048); 
            }
            //<<debug
            for (int j = 0; j < length(seg_clips); j++)
            {
                int bam_start_x = seg_clips_src1[j].first;
                int bam_start_y = seg_clips_src2[j].first;
                int bam_strand = get_cord_strand(seg_str_cord);
                setClippedPositions(row3, row4, seg_clips[j].first, seg_clips[j].second);
                insertNewBamRecord(bam_records, row3, row4, g_id, bam_start_x, bam_start_y, bam_strand, -1, 1, 2048); 
            }
        }
    }
    return 0;
}
int align_gaps (String<BamAlignmentRecordLink> & bam_records,
                GapRecords & gaps,
                StringSet<String<Dna5> >& genomes,
                String<Dna5> & read, 
                String<Dna5> & comrev_read,
                Score<int> & score_scheme,
                AlignGapParms & align_gap_parms)
{
    GapRecordHolder gap(gaps);
    while (!gap.atEnd()) 
    { 
        align_gap(gap, bam_records, genomes, read, 
                comrev_read, score_scheme, align_gap_parms);
        gap.next();
    }
    return 0;
}
/*
 * debug utility
 */
void printGaps(String<std::pair<uint64_t, uint64_t> > & gaps)
{
    for (int i = 0; i < length(gaps); i++)
    {
        std::cout << "[]::gaps " 
                  << get_cord_id(gaps[i].first) << " " 
                  << get_cord_x(gaps[i].first) << " " 
                  << get_cord_y(gaps[i].first) << " " 
                  << get_cord_x(gaps[i].second) << " " 
                  << get_cord_x(gaps[i].second) << " " 
                  << get_cord_y(gaps[i].second) << "\n"; 
    }
}

/*========================================================
=            Main functions of aligning cords            =
=========================================================*/
int const flag_clip_unset = 1 << 32;
int const flag_clip_head = 1;
int const flag_clip_tail = 2;
void set_clip_head_flag (int &flag) {flag |= flag_clip_head;}
void set_clip_tail_flag (int &flag) {flag |= flag_clip_tail;}
int is_clip_head_set (int flag){return flag & flag_clip_head;}
int is_clip_tail_set (int flag){return flag & flag_clip_tail;}

void setLeftClosed(int & f_region)
{
    f_region = -1;
}
void setRghtClosedd(int & f_region)
{
    f_region = 1;
}
void setAllClosed(int & f_region)
{
    f_region = 0; 
}
bool isLeftClosed(int const & f_region)
{
    return f_region == -1;
}
bool isRghtClosed(int const & f_region)
{
    return f_region == 1;
}
bool isAllClosed(int const & f_region)
{
    return f_region == 0;
}
struct AlignCords
{
    String<uint64_t> cords_str;
    String<uint64_t> cords_end;
    String<uint64_t> bands;
    String<uint64_t> bands_lower;
    String<uint64_t> bands_upper;
    String<int> status;
    //if head of alignment of cords_str[i] cords_end[i] need to be clipped
    int if2ClipHead(int i);
    int if2ClipTail(int i);
    //set head of alignment of cords_str[i] cords_end[i] to be clipped
    void set2ClipHead(int i);
    void set2ClipTail(int i);
    int createAlignCords(StringSet<String<Dna5> >& genomes,
                String<Dna5> & read, 
                String<Dna5> & comrev_read,
                String<uint64_t> & cords_str_map,
                String<uint64_t> & cords_end_map,
                int band_lower,
                int band_upper,
                float thd_err_rate,
                int thd_min_abort_anchor);

};
int AlignCords::if2ClipHead(int i)
{
    return is_clip_head_set(status[i]);
}
int AlignCords::if2ClipTail(int i)
{
    return is_clip_tail_set(status[i]);
}
//set head of alignment of cords_str[i] cords_end[i] to be clipped
void AlignCords::set2ClipHead(int i)
{
    set_clip_head_flag(status[i]);
}
void AlignCords::set2ClipTail(int i)
{
    set_clip_tail_flag(status[i]);
}

/**
 * To customize cords to be aligned 
 */
int _initAlignCords(StringSet<String<Dna5> >& genomes,
                    String<Dna5> & read, 
                    String<Dna5> & comrev_read,
                    String<uint64_t> & cords_str_map,
                    String<uint64_t> & cords_end_map,
                    String<uint64_t> & cords_str,
                    String<uint64_t> & cords_end,
                    float thd_err_rate,
                    int thd_min_abort_anchor)
{
    int thd_drop_gap = 20; //drop all gap tiles of record if num of gap blocks > 
    uint64_t recd;
    int recd_str, recd_end; 
    if (empty(cords_str_map))
    {
        return 0;
    }
    appendValue(cords_str, cords_str_map[0]);
    appendValue(cords_end, cords_end_map[0]); 
    for (int i = 1; i < length(cords_str_map); i = recd_end)
    {
        recd = get_cord_recd (cords_str_map[i]);
        recd_str = i;
        recd_end = i + 1; 
        int n_bk = 0; //count blocks within recd
        for (; recd_end < length(cords_str_map); recd_end++)
        {
            if (_DefaultCord.isBlockEnd(cords_str_map[recd_end]))
            {
                ++n_bk;
            }
            if (get_cord_recd(cords_str_map[recd_end]) != recd)
            {
                break;
            }
        }
        int f_dg = (n_bk > thd_drop_gap) ? 1 : 0;
        for (int j = recd_str; j < recd_end; j++)
        {
            uint64_t new_cord_str = cords_str_map[j];
            uint64_t new_cord_end = cords_end_map[j];
            uint64_t cordy = get_cord_y(new_cord_str);
            uint64_t cordx = get_cord_x(new_cord_str);
            int g__id = get_cord_id(new_cord_str);
            int f_drop = 0;
            if (cordy > length(read) - 1 || cordx > length(genomes[g__id]) - 1)
            {
                f_drop = 1;
            }
            else if (!_DefaultCord.isBlockEnd(back(cords_str)) && 
                     !get_cord_strand(back(cords_str) ^ new_cord_str))
            {
                int64_t cordx1 = get_cord_x(back(cords_str));
                int64_t cordx2 = get_cord_x(new_cord_str);
                int64_t cordy1 = get_cord_y(back(cords_str));
                int64_t cordy2 = get_cord_y(new_cord_str);
                if (cordx1 > cordx2 || cordy1 > cordy2) 
                {
                    int64_t danchor =  std::abs(cordx1 - cordy1) -
                                       std::abs(cordx2 - cordy2);
                    if (std::abs(danchor) > 
                        thd_err_rate * std::abs(cordx1 - cordx2) &&
                        std::abs(danchor) > thd_min_abort_anchor) 
                    {
                        set_cord_block_end(back(cords_str));
                        set_cord_block_end(back(cords_end));
                    }
                    else
                    {
                        if (_DefaultCord.isBlockEnd(new_cord_str))
                        {
                            set_cord_block_end(back(cords_str));
                            set_cord_block_end(back(cords_end));
                        }
                        f_drop = 1;
                    }
                }
                else if (get_cord_y(back(cords_end)) < get_cord_y(new_cord_str) ||
                         get_cord_x(back(cords_end)) < get_cord_x(new_cord_str))
                {
                    int64_t dy = get_cord_y(new_cord_str) - get_cord_y(back(cords_end));
                    int64_t dx = get_cord_x(new_cord_str) - get_cord_y(back(cords_end));
                    int64_t thd_cscs_shift = 24;
                    int64_t thd_cscs_same_anchor = 50;
                    if (std::abs(dy - dx) < thd_cscs_same_anchor)
                    {
                        int64_t d = std::max(dy, dx) + thd_cscs_shift;
                        shift_cord(back(cords_end), d, d);
                    }
                    else
                    {
                        int64_t d = std::min(dy, dx) + thd_cscs_shift;
                        shift_cord(back(cords_end), d, d);
                        shift_cord(new_cord_str, -d, -d);
                    }
                }
            }
            if (!f_drop && (is_cord_main(new_cord_str) || !f_dg))
            {
                //if (!mergeAlignCords(back(cords_str), back(cords_end), new_cord_str,
                //      new_cord_end));
                {
                    appendValue (cords_str, new_cord_str);
                    appendValue (cords_end, new_cord_end);
                }
            }
            if (_DefaultCord.isBlockEnd(new_cord_str))
            {
                set_cord_block_end(back(cords_str));
                set_cord_block_end(back(cords_end));
            }
        }
    }
    //print_cords(cords_str, "trmc1");
    return 0;
}
int AlignCords::createAlignCords(StringSet<String<Dna5> >& genomes,
                                 String<Dna5> & read, 
                                 String<Dna5> & comrev_read,
                                 String<uint64_t> & cords_str_map,
                                 String<uint64_t> & cords_end_map,
                                 int band_lower,
                                 int band_upper,
                                 float thd_err_rate,
                                 int thd_min_abort_anchor)
{
    //:Filter
    _initAlignCords(genomes, read, comrev_read, cords_str_map, cords_end_map, 
        cords_str, cords_end, thd_err_rate, thd_min_abort_anchor);
    resize(bands_lower, length(cords_str));
    resize(bands_upper, length(cords_end));
    for (unsigned i = 0; i < length(cords_str); i++)
    {
        bands_upper[i] = band_lower;
        bands_lower[i] = band_upper;
    }
    //:Merge cords and bands
    MergeCordsBandsParm mcb_parms;
    mergeCordsBands(cords_str, cords_end, bands_lower, bands_upper, mcb_parms);
    //:Specific cord_str and cord_end for align
    //_createAlignCords();
    return 0;
}
int check_cord_1_(uint64_t cord, unsigned lx, unsigned ly)
{
    return get_cord_x(cord) >= lx || get_cord_y(cord) >= ly;
}
// return 1 for same strand cord if the pre-cord is larger 
int check_cord_2_(uint64_t cord1, uint64_t cord2)
{
    if (cord1 == cord2 || get_cord_x(cord1) == get_cord_x(cord2) || get_cord_y(cord1) == get_cord_y(cord2))
    {
        return 1;
    }
    if (!_DefaultCord.isBlockEnd(cord1) && !get_cord_strand(cord1 ^ cord2))
    {
        if (get_cord_y(cord1) > get_cord_y(cord2) || 
            get_cord_x(cord1) > get_cord_x(cord2))
        {
            return 1;
        }
    }
    return 0;
}
int check_align_(Row<Align<String<Dna5>, ArrayGaps> >::Type & row1,
                 Row<Align<String<Dna5>, ArrayGaps> >::Type & row2,
                 int score_align,
                 int cord_type_flag,
                 int min_window,
                 int min_score
                 //int g_end //view coordinates
                )
{
    int b1 = beginPosition(row1), b2 = beginPosition(row2);
    int e1 = endPosition(row1), e2 = endPosition(row2);
    float score_density = (float)score_align / std::max(e2 - b2, e1 - b1);
    float score_window = std::max((float) (e1 - b1) / length(getValue(row1._source)), (float)(e2 - b2) / length(getValue(row2._source)));
    if (cord_type_flag == 0) //normal
    {
        if (score_density < s_score_density_thd || score_window < s_score_window_thd) 
        {
            return 1; 
        }
    }
    if (cord_type_flag == 1 || cord_type_flag == -1) //tail clip
    {
        if (e1 - b1 < min_window)
        {
            return 1;
        }
        int sc = 0;
        for (int i = 0; i < min_window; i++)
        {
            if (row1[i] == row2[i]) 
            {
                sc++;
            }
        }
        for (int i = 0; i < e1 - b1 - min_window; i++)
        {
            if (sc > min_score)
            {
                return 0;
            }
            row1[i]==row2[i] ? sc-- : sc;
            row1[i + min_window]==row2[i + min_window] ? sc++ : sc;
        }
        return 1;
    }
    //dout << "ca1" << "\n";
    return 0;
}
/**
 * Main function to align cords and generate bam records.
 */
int alignCords (StringSet<String<Dna5> >& genomes,
                String<Dna5> & read, 
                String<Dna5> & comrev_read,
                String<uint64_t> & cords_str_map, //cords of apx map
                String<uint64_t> & cords_end_map,
                String<BamAlignmentRecordLink> & bam_records,
                int block_size,
                int band) 
{
    double alg_time = sysTime();
    if (length(cords_str_map) < 2) // empty cords
    {
        return 0;
    }
    GapRecords gaps;
    String<TRow> rstr;
    resize(rstr, 8);
    int thd_merge_gap = block_size / 3; //adjacent gaps will be merge if < 
    int thd_max_dshift = block_size * 3; //New regions to be aligned will be extended at two ends to keep overlappad region to be clippped.
    int thd_min_window = 50;
    int thd_min_score = 40;
    float thd_err_rate = 0.25;
    int thd_min_abort_anchor = 20;
    int head_end = block_size >> 2, tail_start = block_size - (block_size >> 2);
    int ri = 0, ri_pre = length(rstr) - 2;
    int flag = 0, flag_pre = 0;
    int f_gap_merge = 0;
    int64_t thd_d_anchor = 50;
    clear (bam_records);

    uint64_t cord_str;
    uint64_t cord_end;
    uint64_t pre_cord_str;
    uint64_t pre_cord_end;
    uint64_t gap_str_cord;
    uint64_t gap_end_cord;
    AlignCords align_cords; 
    double t5 = sysTime();
    align_cords.createAlignCords(genomes, read, comrev_read, cords_str_map, cords_end_map, 
        band, band, thd_err_rate, thd_min_abort_anchor);
    t5 = sysTime() - t5;
    String<uint64_t> & cords_str = align_cords.cords_str;
    String<uint64_t> & cords_end = align_cords.cords_end;
    String<uint64_t> & bands = align_cords.bands;
    double st2 = 0;
    double st3 = sysTime();
    double st4 = 0;
    for (int i = 1; i < (int)length(cords_str); i++)
    {
        int f_clip_head = 1;
        int f_clip_tail = 1;
        int check_flag = 0;
        int g_id = get_cord_id(cords_str[i]);
        uint64_t strand = get_cord_strand(cords_str[i]);
        cord_str = cords_str[i];
        //cord_end = shift_cord(cord_str, block_size, block_size);
        cord_end = cords_end[i];
        if (_DefaultCord.isBlockEnd(cords_str[i - 1])) 
        {
            flag = flag_pre = 0;
            pre_cord_str = pre_cord_end = emptyCord;
            gap_str_cord = gap_end_cord = emptyCord;
            f_gap_merge = 0;
            int dy = std::min(thd_max_dshift, (int)get_cord_y(cords_str[i]));
            int dx = dy;
            cord_str = shift_cord (cords_str[i], -dx, -dy);
            check_flag = -1;
        }
        else if (get_cord_strand (cords_str[i] ^ cords_str[i - 1]))
        {
            int dx = std::min ((int)get_cord_x (cords_str[i] - cords_str[i - 1]), thd_max_dshift);
            int dy = std::min ((int)get_cord_y(cords_str[i]), dx);
            cord_str = _DefaultCord.shift(cords_str[i], -dx, -dy);
            f_clip_head = 0;
            f_clip_tail = 1;
            check_flag = -1;
        }
        if (_DefaultCord.isBlockEnd(cords_str[i]) || i == length(cords_str) - 1)
        {
            int dy = length(read) - get_cord_y(cord_end) - 1;
            dy = std::min(dy, thd_max_dshift);
            int dx = dy;
            cord_end = shift_cord (cord_end, dx, dy);
            check_flag = 1;
            //dout << "shift1" << get_cord_y(cord_str) << get_cord_y(cord_end) << dy << "\n";
        }
        else if (get_cord_strand(cords_str[i] ^ cords_str[i + 1]))
        {
            int dx = std::min((int)get_cord_x (cords_str[i + 1] - cords_str[i]), thd_max_dshift);
            int dy = std::min(dx, int(length(read) - get_cord_y(cords_str[i])));
            cord_end = _DefaultCord.shift(cord_end, dx, dy);
            f_clip_head = 1;
            f_clip_tail = 0;
            check_flag = 1;
        }
        if (i > 1 && 
            std::abs(int64_t(get_cord_x(cords_str[i]) - get_cord_x(cords_str[i - 1]) - 
            get_cord_y(cords_str[i]) + get_cord_y(cords_str[i - 1]))) > thd_d_anchor &&
            !get_cord_strand(cords_str[i] ^ cords_str[i - 1]))
        { //ins, del right cord
            cord_str = shift_cord(cord_str, -50, -50);
            //print_cord(cord_str, "inshit1");
            f_clip_head = 0;
            check_flag = -1; 
        }
        if (i + 1 < length(cords_str) && std::abs(int64_t(get_cord_x(cords_str[i]) - get_cord_x(cords_str[i + 1]) - 
            get_cord_y(cords_str[i]) + get_cord_y(cords_str[i + 1]))) > thd_d_anchor &&
            !get_cord_strand(cords_str[i] ^ cords_str[i + 1]))
        {//ins, del left cord
            cord_end = shift_cord(cord_end, 50, 50);
            f_clip_tail = 0;
            //print_cord(cord_end, "inshit2");
            check_flag = 1;
        }
        if (check_cord_1_(cord_str, length(genomes[g_id]), length(read)) ||
            check_cord_1_(cord_end, length(genomes[g_id]), length(read)) ||
            check_cord_2_(cord_str, cord_end))  
        {
            if (pre_cord_str != emptyCord)
            {
                gap_str_cord = shift_cord(pre_cord_str,
                                          beginPosition(rstr[ri_pre]),
                                          beginPosition(rstr[ri_pre + 1]));
            }
            continue;
        }
        double st1 =  sysTime();
        int score_align = align_cord (rstr[ri], rstr[ri + 1], genomes[g_id], 
                                      read, comrev_read, 
                                      cord_str, cord_end, band, band);
        st2 += sysTime() - st1;

        int f1 = 0, f2 = 0, f3 = 0;
        if (f_clip_head)
        {
            f1 = clip_head_ (rstr[ri], rstr[ri + 1], head_end);
        }
        if (f_clip_tail)
        {
            f2 = clip_tail_ (rstr[ri], rstr[ri + 1], tail_start);
        }
        if (f_clip_head && f_clip_tail)
        {
            f3 = check_align_(rstr[ri], rstr[ri + 1], score_align, check_flag, thd_min_window, thd_min_score);
        }
        flag = f1 | f2 | f3;
        if (flag)
        {
            continue; //alignment quality check, drop poorly aligned 
        }
        uint64_t cord_str_before_merge = cord_str;
            uint64_t tmp1 = pre_cord_str, tmp2= cord_str;
        //<<feature
        String<int> flags;
        String<int> f_merges; //$flags of if call merger_align_
        String<uint64_t> split_cords_str;
        String<uint64_t> split_cords_end;
        uint64_t thd_joint_view_size = block_size + 500000;
        uint64_t thd_split = 3 * thd_joint_view_size;
        if (clippedEndPosition(rstr[ri]) - clippedBeginPosition(rstr[ri]) > thd_split)
        {
            //$This is to split the block of alignment into two sub-blocks
            //$init Global vars
            uint64_t original_clipped_end = clippedEndPosition(rstr[ri]);
            uint64_t split_cord_str = cord_str, split_cord_end = cord_end;
            _DefaultHit.unsetBlockEnd(split_cord_str);
            _DefaultHit.unsetBlockEnd(split_cord_end);

            //$init the first block:head
            appendValue(flags, flag);
            appendValue(f_merges, 1);
            appendValue(split_cords_str, cord_str);
            appendValue(split_cords_end, cord_end);
            setClippedEndPositions(rstr[ri], rstr[ri + 1], 
                clippedBeginPosition(rstr[ri]) + thd_joint_view_size);

            //$init the second block:middle
            
            appendValue(flags, 0);//defined as always successfully merged
            appendValue(f_merges, 0);
            int ri_next = (ri + 2) % length(rstr);
            copyRow(rstr[ri_next], rstr[ri]);
            copyRow(rstr[ri_next + 1], rstr[ri + 1]);
            appendValue(split_cords_str, split_cord_str);
            appendValue(split_cords_end, split_cord_end);
            setClippedBeginPositions(rstr[ri_next], rstr[ri_next + 1], 
                clippedBeginPosition(rstr[ri_next]) + thd_joint_view_size);
            setClippedEndPositions(rstr[ri_next], rstr[ri_next + 1], 
                original_clipped_end - thd_joint_view_size);
            
            //$init the third block:tail 
            appendValue(flags, 0);
            appendValue(f_merges, 0);
            ri_next = (ri + 4) % length(rstr);
            copyRow(rstr[ri_next], rstr[ri]);
            copyRow(rstr[ri_next + 1], rstr[ri + 1]);
            if (_DefaultCord.isBlockEnd(cord_str))
            {
                _DefaultHit.setBlockEnd(split_cord_str);
                _DefaultHit.setBlockEnd(split_cord_end);
            }
            appendValue(split_cords_str, split_cord_str);
            appendValue(split_cords_end, split_cord_end);
            setClippedBeginPositions(rstr[ri_next], rstr[ri_next + 1], 
                original_clipped_end - thd_joint_view_size);
            setClippedEndPositions(rstr[ri_next], rstr[ri_next + 1], original_clipped_end);
            
        }
        else
        {   
            //$Nothing changed when it's not splitted 
            appendValue(flags, flag);
            appendValue(f_merges, 1);
            appendValue(split_cords_str, cord_str);
            appendValue(split_cords_end, cord_end);
        }
        //>>feature 
        for (unsigned j = 0; j < length(split_cords_str); j++)
        {
            cord_str = split_cords_str[j];
            cord_end = split_cords_end[j];
            if (_DefaultCord.isBlockEnd(pre_cord_str))
            {
                //clip_segs(rstr[ri], rstr[ri + 1], cord_str, _gap_parm, -1);
                uint64_t bam_flag = i == 1 ? 0 : 2048;
                insertNewBamRecord(bam_records, g_id,
                                get_cord_x(cord_str) + beginPosition(rstr[ri]),
                                get_cord_y(cord_str) + beginPosition(rstr[ri + 1]),
                                get_cord_strand(cords_str[i]),
                                -1, 1, bam_flag);
                pre_cord_str = cord_str;
                pre_cord_end = cord_end;
                flag = 0;
                flag_pre = 0;
                ri_pre = ri;
                ri = (ri + 2) % length(rstr);
                //std::swap (ri, ri_pre); 
                continue;
            } 
            else if (!f_merges[j])
            {
                flag = flags[j];
            }
            else
            {
                st1 = sysTime();
                flag = merge_align_(rstr[ri_pre], rstr[ri_pre + 1], 
                        rstr[ri], rstr[ri + 1], genomes[g_id], read, comrev_read, 
                        pre_cord_str, cord_str );
                st4 += sysTime() - st1;
            }
            uint64_t bam_start_x = get_cord_x(pre_cord_str) + beginPosition(rstr[ri_pre]);
            uint64_t bam_start_y = get_cord_y(pre_cord_str) +  beginPosition(rstr[ri_pre + 1]);
            uint64_t bam_strand = get_cord_strand(pre_cord_str); 
            if (!flag_pre)
            {
                if (!flag)
                {
                    insertBamRecordCigar(back(bam_records), 
                                rstr[ri_pre], rstr[ri_pre + 1]);                
                    f_gap_merge = 0;
                }
                else if (flag & 1)
                {
                    gap_str_cord = shift_cord(pre_cord_str,
                                              beginPosition(rstr[ri_pre]),
                                              beginPosition(rstr[ri_pre + 1]));
                }
                else if (flag & 2)
                {
                    //clip_segs(rstr[ri_pre], rstr[ri_pre + 1], 
                    //          pre_cord_str, _gap_parm, 1); 
                    insertBamRecordCigar(back(bam_records), 
                                         rstr[ri_pre], 
                                         rstr[ri_pre + 1]);            
                }
            }
            else if (flag_pre & 1)
            {
                if (!flag)
                {
                    gap_end_cord = shift_cord(pre_cord_str, 
                                              endPosition(rstr[ri_pre]), 
                                              endPosition(rstr[ri_pre + 1]));
                    if (get_cord_x(gap_str_cord) < get_cord_x(gap_end_cord) &&
                        get_cord_y(gap_str_cord) < get_cord_y(gap_end_cord))
                    {
                        if(insertGaps(gaps, gap_str_cord, gap_end_cord,
                                      length(bam_records) - 1,
                                      thd_merge_gap, f_gap_merge))
                        {
                            bam_start_x = get_cord_x(cord_str) + beginPosition(rstr[ri]);
                            bam_start_y = get_cord_y(cord_str) + beginPosition(rstr[ri + 1]);
                            bam_strand = get_cord_strand(cord_str);
                            insertNewBamRecord(bam_records, g_id, bam_start_x, bam_start_y, 
                                bam_strand, -1, 1, 2048); 
                        }      
                        f_gap_merge = 1;
                    }
                }
                else if (flag & 1)
                {
                   //NONE 
                }
                else if (flag & 2)
                {
                    //clip_segs(rstr[ri_pre], rstr[ri_pre + 1], 
                    //          pre_cord_str, _gap_parm, 1); 
                    gap_end_cord = shift_cord(pre_cord_str, 
                                              endPosition(rstr[ri_pre]), 
                                              endPosition(rstr[ri_pre + 1]));
                    if (get_cord_x(gap_str_cord) < get_cord_x(gap_end_cord) &&
                        get_cord_y(gap_str_cord) < get_cord_y(gap_end_cord))
                    {
                        if(insertGaps(gaps, gap_str_cord, gap_end_cord,
                                      length(bam_records) - 1,
                                      thd_merge_gap, f_gap_merge))
                        {

                            bam_start_x = get_cord_x(cord_str) + beginPosition(rstr[ri]);
                            bam_start_y = get_cord_y(cord_str) + beginPosition(rstr[ri + 1]);
                            bam_strand = get_cord_strand(cord_str);   
                            insertNewBamRecord(bam_records, g_id, bam_start_x, bam_start_y, bam_strand, -1, 1, 2048); 
                        }               
                        f_gap_merge = 1;
                    }
                }
            }
            else if (flag_pre & 2) //diff strands
            {
                if (!flag)
                {
                    //clip_segs(rstr[ri_pre], rstr[ri_pre + 1], 
                    //          pre_cord_str, _gap_parm, -1);
                    bam_start_x = get_cord_x(pre_cord_str) + 
                                  beginPosition(rstr[ri_pre]);
                    bam_start_y = get_cord_y(pre_cord_str) + 
                                  beginPosition(rstr[ri_pre + 1]);
                    insertNewBamRecord(bam_records, 
                                       rstr[ri_pre], 
                                       rstr[ri_pre + 1],
                                       g_id, bam_start_x, bam_start_y, bam_strand,
                                       -1, 1, 2048);
                    f_gap_merge = 0;                     

                }
                else if (flag & 1)
                {
                    //clip_segs(rstr[ri_pre], rstr[ri_pre + 1], 
                    //          pre_cord_str, _gap_parm, -1);     
                    gap_str_cord = shift_cord(pre_cord_str,
                                              beginPosition(rstr[ri_pre]),
                                              beginPosition(rstr[ri_pre + 1]));
                    bam_start_x = get_cord_x(pre_cord_str) + 
                                  beginPosition(rstr[ri_pre]);
                    bam_start_y = get_cord_y(pre_cord_str) +
                                  beginPosition(rstr[ri_pre + 1]);
                    insertNewBamRecord(bam_records, 
                                       g_id, bam_start_x, bam_start_y, bam_strand,
                                       -1, 1, 2048);
                    f_gap_merge = 0;
                }
                else if (flag & 2)
                {
                    //clip_segs(rstr[ri_pre], rstr[ri_pre + 1], 
                    //          pre_cord_str, _gap_parm, 0); 
                    bam_start_x = get_cord_x(pre_cord_str) +
                                  beginPosition(rstr[ri_pre]);
                    bam_start_y = get_cord_y(pre_cord_str) +
                                  beginPosition(rstr[ri_pre + 1]);
                    insertNewBamRecord(bam_records, 
                                       rstr[ri_pre], 
                                       rstr[ri_pre + 1],
                                       g_id, bam_start_x, bam_start_y, bam_strand,
                                       -1, 1, 2048);  
                    f_gap_merge = 0;
                }
            }
            //addition process for the last cord of block 
            if (_DefaultCord.isBlockEnd(cords_str[i]))
            {
                if (!flag)
                {
                    //todo::clipping last 30 bases.
                    //clip_segs(rstr[ri], rstr[ri + 1], 
                    //          cord_str, _gap_parm, 1); 
                    insertBamRecordCigar(back(bam_records), rstr[ri], rstr[ri + 1]);                 
                }
                else if (flag & 1)
                {
                    //For similicity just clip and insert new bam
                    //but todo::better to add gaps.
                    //clip_segs(rstr[ri], rstr[ri + 1], 
                    //          cord_str, _gap_parm, 0); 
                    bam_start_x = get_cord_x(cord_str) +
                                  beginPosition(rstr[ri]);
                    bam_start_y = get_cord_y(cord_str) +
                                  beginPosition(rstr[ri + 1]);
                    bam_strand = get_cord_strand (cord_str);
                    insertNewBamRecord(bam_records, 
                                       rstr[ri], 
                                       rstr[ri + 1],
                                       g_id, bam_start_x, bam_start_y, bam_strand,
                                       -1, 1, 2048);
                            //dout << "ib14" << bam_start_y << "\n";

                }
                else if (flag & 2)
                {
                    //clip_segs(rstr[ri], rstr[ri + 1], 
                    //          cord_str, _gap_parm, 0); 
                    bam_start_x = get_cord_x(cord_str) +
                                  beginPosition(rstr[ri]);
                    bam_start_y = get_cord_y(cord_str) +
                                  beginPosition(rstr[ri + 1]);
                    bam_strand = get_cord_strand(cord_str);
                    insertNewBamRecord(bam_records, 
                                       rstr[ri], 
                                       rstr[ri + 1],
                                       g_id, bam_start_x, bam_start_y, bam_strand,
                                       -1, 1, 2048);               
                                                    //dout << "ib12" << bam_start_y << "\n";

                }
                flag_pre = 0;
                flag = 0;
                pre_cord_str = pre_cord_end = emptyCord;
                cord_str = cord_end = emptyCord;
                f_gap_merge = 0;
            }
            else
            {
                flag_pre = flag;
                flag = 0;
                pre_cord_str = cord_str;
                pre_cord_end = cord_end;
                ri_pre = ri;
                ri = (ri + 2) % length(rstr);
                //std::swap (ri, ri_pre); //swap the current and pre row id in the aligner.
            }
        }
    }
    //printGaps(gaps.c_pairs);
    Score<int> score_scheme;
    /*
    int thd_clip_score = 80;
    int thd_reject_score = 130;
    int thd_accept_score = 140;
    int thd_accept_density = 16;
    */
    align_gaps(bam_records, gaps, genomes, read, comrev_read, score_scheme, _gap_parm); 
    clipBamRecordLinksHeadTail(bam_records, length(read), _gap_parm);
    //printCigarSrcLen(bam_records, "pscr_gggaps1 ");
        
    return 0;
}
/*=====  End of Main functions of aligning cords  ======*/
