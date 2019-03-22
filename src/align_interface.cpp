#include <utility> 
#include "pmpfinder.h"
#include "align_interface.h"
//TODO seqand::setclippedpositoin retrieve source postion that's not efficient
//TODO make holder to rows, so set clip postion can be iteratored.
using namespace seqan;
/**
 * seqan::view coordinate and source coordinate transformation NOTES
 * setClippedBeginPosition(row1, c1)
 * @c1 := current_view_coordinate + clippedBeginPosition(row1)
 */
typedef Align<String<Dna5>, ArrayGaps> TAlign;
typedef Row<TAlign>::Type TRow; 
typedef Iterator<TRow>::Type TRowIterator;

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
    std::cout << "atEnd " << length(holder.c_pairs) << "\n";
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
int const s3 = -1; //gap
float s_score_density_thd = 2; //if < the value alignment of cords will be dropped
float s_score_window_thd = 0.75;
int thd_align_score = 350 /*depends on score_scheme*/;

uint16_t bam_flag_rvcmp = 16;
uint16_t bam_flag_rvcmp_nxt = 32;
uint16_t bam_flag_suppl = 2048;
uint64_t CORD_NULL = _DefaultCord.makeBlockEndVal(~0);
uint64_t emptyCord = CORD_NULL - 1;

GapParm _gap_parm;
Score<int, Simple> _default_scheme_ (s1, s2, s3);

int const flag_clip_unset = 1 << 32;
int const flag_clip_head = 1;
int const flag_clip_tail = 2;
int const flag_clip_both = flag_clip_head | flag_clip_tail;
void set_clip_head_flag (int &flag) {flag |= flag_clip_head;}
void set_clip_tail_flag (int &flag) {flag |= flag_clip_tail;}
int is_clip_head_set (int flag){return flag & flag_clip_head;}
int is_clip_tail_set (int flag){return flag & flag_clip_tail;}


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
GapParm::GapParm ():
        thd_clip_score(80),   //density score: clip when density is < 0.8
        thd_reject_score(130),      //calculate density within the window
        thd_accept_score(140),  //the mini distanse between two clips.
        thd_min_interval(20),
        thd_accept_density(4.5)
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

/*
 * insert cigar to the original cigar 
 */
int insertBamRecordCigar (BamAlignmentRecord & bam_record,
                    Row<Align<String<Dna5>, ArrayGaps> >::Type & row1,
                    Row<Align<String<Dna5>, ArrayGaps> >::Type & row2, 
                    int pos = -1
                   )
{
    //<<debug
    std::cout <<"xxxxxxxxxxxxxxxxinscigar\n";
    std::pair <int,int> lens;
    lens = countCigar (bam_record.cigar);
    std::cout << "cigar_before " << lens.first << " " << lens.second << "\n";
    //>>debug
    if (pos < 0)
    {
        align2cigar(bam_record.cigar, row1, row2);
    }
    else
    {
        if (pos > length(bam_record.cigar ) - 1)
        {
            //<<debug
            lens = countCigar (bam_record.cigar);
            std::cout << "cigar_after " << lens.first << " " << lens.second << "\n"; 
            //>>debug
            return 1;
        }
        String<CigarElement< > > tmp;
        align2cigar(tmp, row1, row2);
        std::cout << "inscigar " << length(tmp) << "\n";
        insertCigar(bam_record.cigar, pos, tmp);
    }
    //<<debug
    lens = countCigar (bam_record.cigar);
    std::cout << "cigar_after " << lens.first << " " << lens.second << "\n";
    //>>debug
    return 0;
}
int  insertNewBamRecord (String<BamAlignmentRecordLink> & bam_records,
                    Row<Align<String<Dna5>, ArrayGaps> >::Type & row1,
                    Row<Align<String<Dna5>, ArrayGaps> >::Type & row2, 
                    int g_id,
                    int g_beginPos,
                    int strand,
                    int pos = -1
                   )
{
    BamAlignmentRecordLink bam_record;
    insertBamRecordCigar(bam_record, row1, row2);
    if (g_id >= 0)
    {
        bam_record.rID = g_id;
    }
    if (g_beginPos >= 0)
    {
        bam_record.beginPos = g_beginPos; 
    }
    if (strand)
    {
        bam_record.flag |= bam_flag_rvcmp | bam_flag_suppl;
    }
    if (pos < 0)
    {
        appendValue(bam_records, bam_record);
    }
    else 
    {
        insert(bam_records, pos, bam_record);
    }
    return 0;
}
int  insertNewBamRecord (String<BamAlignmentRecordLink> & bam_records,
                                int g_id,
                                int g_beginPos,
                                int strand,
                                int pos = -1
                               )
{
    BamAlignmentRecordLink bam_record;
    if (g_id >= 0)
    {
        bam_record.rID = g_id;
    }
    if (g_beginPos >= 0)
    {
        bam_record.beginPos = g_beginPos; 
    }
    if (strand)
    {
        bam_record.flag |= bam_flag_rvcmp | bam_flag_suppl;
    }
    if (pos < 0)
    {
        appendValue(bam_records, bam_record);
    }
    else 
    {
        insert(bam_records, pos, bam_record);
    }
    return 0;
}
/*
 * Cigar from the row1 and row2 are inserted to the cigar from the bam_record
 */
int  insertBamRecord (BamAlignmentRecord & bam_record,
                    Row<Align<String<Dna5>, ArrayGaps> >::Type & row1,
                    Row<Align<String<Dna5>, ArrayGaps> >::Type & row2, 
                    int g_id,
                    int g_beginPos,
                    int pos = -1
                   )
{
    insertBamRecordCigar(bam_record, row1, row2, pos);
    if (g_id >= 0)
    {
        bam_record.rID = g_id;
    }
    if (g_beginPos >= 0)
    {
        bam_record.beginPos = g_beginPos; 
    }
    return 0;
}
/**
 * debug utility
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
void printCigar(String<CigarElement< > > &cigar)
{
    for (int i = 0; i < length(cigar); i++)
    {
        std::cout << cigar[i].count <<cigar[i].operation;
    }
    std::cout << "\n";
}
 int align_block (Row<Align<String<Dna5>, ArrayGaps> >::Type & row1,
                        Row<Align<String<Dna5>, ArrayGaps> >::Type & row2,
                        String<Dna5> & genome,
                        String<Dna5> & read,
                        String<Dna5> & comrevRead,
                        uint64_t strand,
                        uint64_t genomeStart,
                        uint64_t genomeEnd,
                        uint64_t readStart,
                        uint64_t readEnd,
                        int band,
                        Score<int> scheme = _default_scheme_
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
    int score = globalAlignment(row1, row2, scheme, AlignConfig<true, true, true, true>(), -band, band);
    return score;
}

int cord2row_(Row<Align<String<Dna5>, ArrayGaps> >::Type & row1,
                Row<Align<String<Dna5>, ArrayGaps> >::Type & row2,
                String<Dna5> & genome,
                String<Dna5> & read, 
                String<Dna5> & comrevRead,
                uint64_t cord_start,
                uint64_t cord_end
               )
{
    uint64_t genomeStart = get_cord_x(cord_start);
    uint64_t genomeEnd = get_cord_x(cord_end);
    uint64_t readStart = get_cord_y(cord_start);
    uint64_t readEnd = get_cord_y(cord_end);
    uint64_t strand = get_cord_strand (cord_start);   
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
                uint64_t cord_start,
                uint64_t cord_end,
                int band,
                int local_flag = 1,
                Score<int> scheme = _default_scheme_
               )
{
    cord2row_ (row1, row2, genome, read, comrevRead, cord_start, cord_end);
    int score = 0;
    if (!local_flag)
    {
        score = localAlignment (row1, row2, scheme, DynamicGaps());
    }
    else
    {
        score = globalAlignment (row1, row2, scheme, AlignConfig<true, true, true, true>(), -band, band);
    }
    return score;
}
int align_cord (Row<Align<String<Dna5>, ArrayGaps> >::Type & row1,
                Row<Align<String<Dna5>, ArrayGaps> >::Type & row2,
				String<Dna5> & genome,
                String<Dna5> & read, 
                String<Dna5> & comrevRead,
                uint64_t & cord,
                int block_size = window_size,
                int band = window_size / 2,
                int local_flag = 1,
                Score<int> scheme = _default_scheme_
               )
{
    uint64_t cord_end = _DefaultCord.shift(cord, block_size, block_size);
    int score = align_cord (row1, row2, genome, read, comrevRead, cord, cord_end, band, local_flag, _default_scheme_);
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
    int r_flag = 0;
    if (cord1 == emptyCord)
    {
        return 32;
    }
    int thd_cord_overlap = std::max(endPosition(row11), endPosition(row12));
    if (get_cord_strand(cord1 ^ cord2))
    {
        r_flag |= 2;
        return r_flag;
    }
    else if (!_DefaultCord.isCordsOverlap(cord1, cord2, thd_cord_overlap))  //sv: gap or reverse
    {
        r_flag |= 64 | 1;
    }
    if (r_flag)
    {
        return r_flag;
    } 
    int bit = 20, bit2 = 40;
    uint64_t start11 = get_cord_x(cord1);
    uint64_t start21 = get_cord_x(cord2);
    uint64_t start12 = get_cord_y(cord1);
    uint64_t start22 = get_cord_y(cord2);
    int64_t mask = (1ULL << bit) - 1;
	int64_t delta1 = start21 - start11;
	int64_t delta2 = start22 - start12;
	String<int64_t> align1, align2; 
    TRowIterator it1, it2;
    if (endPosition(row11) < beginPosition(row21) + delta1)
    {
        return 1;
    }
    if (endPosition(row12) < beginPosition(row22) + delta2)
    {
        return 1|4096;
    }
    if (beginPosition(row11) > beginPosition(row21) + delta1)
    {
        return 1|128;
    }
    if (beginPosition(row12) > beginPosition(row22) + delta2)
    {
        return 1|256;
    }
    if (endPosition(row11) > endPosition(row21) + delta1)
    {
        return 1|1024;
    }
    if (endPosition(row12) > endPosition(row22) + delta2)
    {
        return 1|2048;
    }
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
        return 8 | 1;
    }
    
    int thd_merge_x = 2, thd_merge_y = 2;
    int flag = 0, start_j = 0;
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
    return 16 | 1;
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
int insertGaps(GapRecords & gaps,
                  uint64_t cord_str, //start cord 
                  uint64_t cord_end, //end cord 
                  int bam_segs_id,
                  int thd_merge_gap
                 )
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
        std::cout << "ig1 " << get_cord_y(back(gaps.c_pairs).second) << " " << get_cord_y(cord1) << "\n";
        if (_DefaultCord.isCordsOverlap(back(gaps.c_pairs).second, cord1, thd_merge_gap) )
        {
        std::cout << "ig2 " << get_cord_y(back(gaps.c_pairs).second) << " " << get_cord_y(cord1) << "\n";
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
                  int thd_merge_gap
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
        if (_DefaultCord.isCordsOverlap(back(gaps.c_pairs).second, cord1, thd_merge_gap) )
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
 void _nextView2Src (Iterator<Row<Align<String<Dna5>, ArrayGaps>>::Type>::Type & it1,
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
 *  The size of the block to be clipped are supposed to > window_size,
 *  since the function clips the middle part of the block excluding 
 *  the head and tail part. Too small block cann't be well clipped. 
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
               int delta = 3
              )
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
        _nextView2Src(it1, it2, src1, src2);
        if (count++ == delta)
        {
            buffer[buf_len] = x; //+ thd_1 * (std::abs(src1 - src1_2 - src2 + src2_2));
            buffer_src1[buf_len] = src1;
            buffer_src2[buf_len] = src2;
            buf_len++;
            count = 1;
        }
        ++it1;
        ++it2;
    }
//TODO::change the score so it will not clip long ins or dels.
    resize (buffer, buf_len);
    int prebp = 0; //pre-breakpoint counter of buffer
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
        std::cout << "clip_gaps_segs::tri " << i << " " 
                    << d_score << " " 
                    << thd_clip_score << " "
                    << d_score_left << " " 
                    << d_score_right << " " 
                    << thd_reject_score << " "
                    << float(d_score_left) / (d_w * delta) << "\n";
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
            float dens_ = (float)(buffer[bp] - buffer[prebp]) / (delta * (bp - prebp + 1));
            std::cout << "clip_gaps_segs::tri <<<<<<<<<<<<<<<,," << i << " " << d_score << " " << dens_ << " " <<thd_clip_score<<"\n";
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
    for (int i = 0; i < length(clip_records); i++)
    {
        std::cout << "[]::clip_records3 " << clip_records[i].first << " " << clip_records[i].second << "\n";
    }
    return 0;
}
int clip_rows_segs (Row<Align<String<Dna5>,ArrayGaps> >::Type & row1,
               Row<Align<String<Dna5>,ArrayGaps> >::Type & row2,
               String<std::pair<int, int> > & clip_records,     //view coordinate
               String<std::pair<int, int> > & clip_records_src1,//source1 coordinate
               String<std::pair<int, int> > & clip_records_src2,//source2 coordinate
               uint64_t gap_start_cord,
               GapParm gap_parm,
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
                    gap_parm.thd_clip_score,
                    gap_parm.thd_reject_score,
                    gap_parm.thd_accept_score,
                    gap_parm.thd_min_interval,
                    gap_parm.thd_accept_density,
                    thd_clip_mini_region,
                    window, //supposed to < 50
                    delta
                   );
}

/*
 * Clip head or tail by calling the clip_gaps_segs. (though not efficient) 
 * @direction -1 head, 1 tail
 */
int clip_segs(Row<Align<String<Dna5>,ArrayGaps> >::Type & row1,
               Row<Align<String<Dna5>,ArrayGaps> >::Type & row2,
               uint64_t gap_start_cord,
               GapParm & gap_parm,
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
                   gap_parm.thd_clip_score,
                   gap_parm.thd_reject_score, 
                   gap_parm.thd_accept_score, 
                   gap_parm.thd_min_interval,
                   gap_parm.thd_accept_density
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

/*
 *!Make sure head and tail in the gap specified by @gap can be well aligned.
 * such that the gap can be merged into the main alignment;
 *@thd_alg_extnd: extension in inversion sub region to be aligned
 * such that it can be clipped precisely at two ends;
 */
int align_gap (GapRecordHolder & gap,
               String<BamAlignmentRecordLink> & bam_records,
               StringSet<String<Dna5> >& genomes,
               String<Dna5> & read, 
               String<Dna5> & comrevRead,
               Score<int> & score_scheme,
               GapParm & gap_parm)
{
    typedef String<std::pair<int, int> > ClipRecords;
    int thd_alg_extnd = 20;
    uint64_t str_cord = gap.getCords().first;
    uint64_t end_cord = gap.getCords().second;
    std::cout << "ag1 " << get_cord_x(str_cord) << " " << get_cord_x(end_cord) << " " 
              << get_cord_y(str_cord) << " " << get_cord_y(end_cord) << "\n";
    if (get_cord_strand(str_cord ^ end_cord))
    {
        return 1;
    }
    ClipRecords clips, clips_src1, clips_src2;
    ClipRecords seg_clips, seg_clips_src1, seg_clips_src2;
    TRow row1, row2 ;
    uint64_t seg_str_cord;
    uint64_t seg_end_cord;
    int g_id = get_cord_id(str_cord);
    int band = std::max(get_cord_x(end_cord - str_cord),
                        get_cord_y(end_cord - str_cord)) / 2;
    int bam_id = gap.getBamSegIdHead();
    int bam_next_id = gap.getBamSegIdTail();
    //WARNING::modify band::too large band
    std::cout << "alg1 " << get_cord_y(str_cord) << " " << get_cord_y(end_cord) << "\n";
    return 0;
    align_cord (row1, row2, genomes[g_id], read, comrevRead, str_cord, end_cord, band);
    printRows(row1, row2, " ag2 ");
    int const view_str = 0;
    int const view_end = clippedEndPosition(row1);
    clip_rows_segs (row1, row2, clips, clips_src1, clips_src2, str_cord, gap_parm);
    if (empty(clips))
    {
        return 1;
    }
    if (length (clips) == 1)
    {
        insertBamRecordCigar(bam_records[bam_id], row1, row2);
        bam_records[bam_id].addNext(bam_next_id);
        std::cout << "ag3 " << " " << get_cord_y(str_cord) << " " << bam_id << " " << bam_next_id << length(bam_records) << "\n";
    }
    else if (length(clips) > 1)
    {
        setClippedPositions(row1, row2, view_str, clips[0].second);
        insertBamRecordCigar(bam_records[bam_id], row1, row2);
        setClippedPositions(row1, row2, back(clips).first, view_end);
        int64_t bam_start = get_cord_x (str_cord) + beginPosition(row1);
        insertBamRecord(bam_records[bam_next_id], row1, row2, g_id, bam_start, 0);
        for (int i = 1; i < length(clips) - 1; i++)
        {
            int bam_start = beginPosition(clips[i].first) + get_cord_x(str_cord);
            int bam_strand = get_cord_strand(str_cord);
            setClippedPositions(row1, row2, clips[i].first, clips[i].second);
            insertNewBamRecord(bam_records, row1, row2, g_id, bam_start, bam_strand); 
        }
        for (int i = 0; i < length(clips) - 1; i++)
        {
            seg_str_cord = set_cord_xy (str_cord, 
                                        clips_src1[i].second, 
                                        clips_src2[i].second);
            seg_end_cord = set_cord_xy (end_cord, 
                                        clips_src1[i + 1].first, 
                                        clips_src2[i + 1].second);
            cmpRevCord (seg_str_cord, seg_end_cord, 
                        seg_str_cord, seg_end_cord, length(read));
            seg_str_cord = shift_cord(seg_str_cord, -thd_alg_extnd, -thd_alg_extnd);
            seg_end_cord = shift_cord(seg_end_cord, thd_alg_extnd, thd_alg_extnd);
            int seg_band = std::max(get_cord_x(seg_end_cord - seg_str_cord),
                                    get_cord_y(seg_end_cord - seg_str_cord)) / 2;
            align_cord (row1, row2, genomes[g_id], 
                        read, comrevRead, seg_str_cord, seg_end_cord, seg_band);
            printRows (row1, row2, " ac4 ");
            clip_rows_segs (row1, row2, seg_clips, seg_clips_src1, seg_clips_src2,
                            seg_str_cord,gap_parm);
            for (int j = 0; j < length(seg_clips); j++)
            {
                int bam_start = seg_clips_src1[j].first;
                int bam_strand = get_cord_strand(seg_str_cord);
                setClippedPositions(row1, row2, seg_clips[i].first, seg_clips[i].second);
                insertNewBamRecord(bam_records, row1, row2, g_id, bam_start, bam_strand); 
                std::cout << "ac4 " << bam_start << " " << seg_clips_src1[j].first << "\n";
            }
        }
    }
    return 0;
}
int align_gaps (String<BamAlignmentRecordLink> & bam_records,
                GapRecords & gaps,
                StringSet<String<Dna5> >& genomes,
                String<Dna5> & read, 
                String<Dna5> & comrevRead,
                Score<int> & score_scheme,
                GapParm & gap_parm)
{
    GapRecordHolder gap(gaps);
    while (!gap.atEnd()) 
    { 
        align_gap(gap, bam_records, genomes, read, 
                comrevRead, score_scheme, gap_parm);
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
    std::cout << "score_density 0 " << score_density << " " << score_window << "\n";
    if (cord_type_flag == 0) //normal
    {
        if (score_density < s_score_density_thd || score_window < s_score_window_thd) 
        {
            std::cout << "score_density 1xxxxxxxxxxxxx " << score_density << " " << "\n";
            return 1; 
        }
    }
    if (cord_type_flag == 1 || cord_type_flag == -1) //tail clip
    {
        std::cout << "check_align_1 \n";
        printRows(row1, row2);
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
    return 0;
}

/*
 *  Align cords and output cigar string.
 *  Steps:
 *  1. Align: Each normal cord will be aligned locally 
 *       First, last cord in the block and cord at two ends of two strands will extend 
 *       the region to be aligned to cover the breakpoints.
 *  2. Trim: only clip consecutive ins or dels at head and tail of the cord in this step) 
 *           and check quality of each normal cord
 *  3. Merge: merge current alignment to its pre alignment.
 *  4. Insert bam and insert gaps:
 *  If cord is normal it (must and only!) clip tail of its pre_align and head of itself. 
 *  If cord is the first in block, it will not clip pre tail but will start a new empty bam record;
 *  If cord is the last in block, it also needs to clip its own tail and add its self to the last bam record;
 *  If cord with its pre are two ends of strands, it will clip pre tail and its own head and start a new empty bam record;
 *  If cord with its pre are two ends of gap, it will insert gap and not clip any.
 *  5. align all gaps.
 */
int align_cords (StringSet<String<Dna5> >& genomes,
                 String<Dna5> & read, 
                 String<Dna5> & comrevRead,
                 String<uint64_t> & cords,
                 String<BamAlignmentRecordLink> & bam_records,
                 int p,
                 int block_size,
                 int band
                ) 
{
    GapRecords gaps;
    BamAlignmentRecordLink emptyBamRecord;
    String<TRow> rstr;
    String<uint64_t> cords_buffer;
    int head_end = block_size >> 2, tail_start = block_size - (block_size >> 2);
    int ri = 0, ri_pre = 2; //cliped segment and row id
    int g_id = -1, g_beginPos = 0, strand = 0, flag = 0, flag_pre = 0;
    int thd_merge_gap = block_size / 3; // two adjacent gaps will be merged to one if < it
    int flag2 = 1;
    int d_overlap_x;
    int d_overlap_y;
    int thd_max_dshift = block_size * 3; //For head and tail cords of a new segs, their size will be increased based on the window_size, so that there can be enough overlap region to be clippped.
    int thd_min_window = 50;
    int thd_min_score = 40;
    float thd_ddx = (float)band / block_size / 8; //TODO::8 needs tunning
    clear (bam_records);
    resize(cords_buffer, 2);
    resize(rstr, 6);
    double t1, t2 = 0, t3 = sysTime();
    if (length(cords) < 2) // cords is empty
    {
        return 0;
    }
    uint64_t cord_start;
    uint64_t cord_end;
    uint64_t pre_cord_start;
    uint64_t pre_cord_end;
    uint64_t gap_str_cord;
    uint64_t gap_end_cord;
    uint64_t bam_start;
    uint64_t bam_strand;
    int check_flag = 0;
    TRow rt1, rt2;
    for (int i = 1; i < (int)length(cords); i++)
    {
        flag_pre = flag;
        flag = 0;
        check_flag = 0;
        g_id = get_cord_id(cords[i]);
        g_beginPos = get_cord_x (cords[i]);
        strand = get_cord_strand(cords[i]);
        cord_start = cords[i];
        cord_end = _DefaultCord.shift(cord_start, block_size, block_size);
        int ri_gap_h = 4;
        if (_DefaultCord.isBlockEnd(cords[i - 1])) 
        {
            int dy = std::min(thd_max_dshift, (int)get_cord_y(cords[i]));
            int dx = dy + thd_ddx * dy;
            pre_cord_start = emptyCord; //merge_align_ return 0 if pre_cord_start is empty
            cord_start = _DefaultCord.shift (cords[i], -dx, -dy);
            check_flag = -1;
            std::cout << "ac2 " << get_cord_y(cord_start) << "\n";
        }
        else if (get_cord_strand (cords[i] ^ cords[i - 1]))
        {
            int dx = std::min ((int)get_cord_x (cords[i] - cords[i - 1]), thd_max_dshift);
            int dy = std::min (dx, (int)get_cord_y(cords[i - 1]));
            cord_start = _DefaultCord.shift(cords[i], -dx, -dy);
            check_flag = -1;
        }
        if (_DefaultCord.isBlockEnd(cords[i]))
        {
            int dy = std::min ((int)length(read) - (int)get_cord_y(cord_end), thd_max_dshift);
            int dx = dy + thd_ddx * dy;
            cord_end = _DefaultCord.shift (cord_end, dx, dy);
            check_flag = 1;
        }
        else if (get_cord_strand(cords[i] ^ cords[i + 1]))
        {
            int dx = std::min((int)get_cord_x (cords[i + 1] - cords[i]), thd_max_dshift);
            int dy = std::min(dx, int(length(read) - get_cord_y(cords[i])));
            cord_end = _DefaultCord.shift(cord_end, dx, dy);
            check_flag = 1;
        }
        int score_align = align_cord (rstr[ri], rstr[ri + 1], genomes[g_id], 
                                      read, comrevRead, 
                                      cord_start, cord_end, band);
        std::cout << "ac1_r " << i + 1 << " " 
        << get_cord_y(cord_start) << " " << get_cord_x(cord_start) <<  "\n";
        std::cout << rstr[ri] << "|" << rstr[ri + 1] << "\n";
        printRows(rstr[ri], rstr[ri + 1], "ac2_2");
        flag = clip_head_ (rstr[ri], rstr[ri + 1], head_end)
             | clip_tail_ (rstr[ri], rstr[ri + 1], tail_start)
             | check_align_(rstr[ri], rstr[ri + 1], score_align, check_flag, thd_min_window, thd_min_score);
        //TODO!!check_align for first or last cord are differnet from middle cords
        if (flag){
            continue; //alignment quality check, drop poorly aligned 
        }
        if (_DefaultCord.isBlockEnd(pre_cord_start))
        {
            clip_segs(rstr[ri], rstr[ri + 1], cord_start, _gap_parm, -1);
            insertNewBamRecord(bam_records, g_id,
                            get_cord_x(cord_start) + beginPosition(rstr[ri]),
                            get_cord_strand(cords[i]));
            pre_cord_start = cord_start;
            pre_cord_end = cord_end;
            flag = 0;
            std::swap (ri, ri_pre); 
            continue;
        } 
        else
        {
            flag = merge_align_(rstr[ri_pre], rstr[ri_pre + 1], 
                    rstr[ri], rstr[ri + 1], pre_cord_start, cord_start );
        }
        std::cout << "flag " << i + 1 << " " << flag_pre << " " << flag << "\n";
        bam_start = get_cord_x(pre_cord_start) + beginPosition(rstr[ri_pre]);
        bam_strand = get_cord_strand(pre_cord_start); 
        if (!flag_pre)
        {
            if (!flag)
            {
                insertBamRecordCigar(back(bam_records), 
                            rstr[ri_pre], rstr[ri_pre + 1]);                
            }
            else if (flag & 1)
            {
                //detach_copy_row(rstr[ri_gap_h], rstr[ri_pre]);
                //detach_copy_row(rstr[ri_gap_h + 1], rstr[ri_pre + 1]);
                gap_str_cord = shift_cord(pre_cord_start,
                                          beginPosition(rstr[ri_pre]),
                                          beginPosition(rstr[ri_pre + 1]));
                std::cout << "ac12 " << get_cord_y(gap_str_cord) << "\n";
            }
            else if (flag & 2)
            {
                clip_segs(rstr[ri_pre], rstr[ri_pre + 1], 
                          pre_cord_start, _gap_parm, 1); 
                insertBamRecordCigar(back(bam_records), 
                                     rstr[ri_pre], 
                                     rstr[ri_pre + 1]);            
            }
        }
        else if (flag_pre & 1)
        {
            if (!flag)
            {
                gap_end_cord = shift_cord(pre_cord_start, 
                                          endPosition(rstr[ri_pre]), 
                                          endPosition(rstr[ri_pre + 1]));
                if(insertGaps(gaps, gap_str_cord, gap_end_cord,
                           //rstr[ri_gap_h], 
                           //rstr[ri_gap_h + 1],
                           //rstr[ri_pre],
                           //rstr[ri_pre + 1],
                           length(bam_records) - 1,
                           thd_merge_gap))
                {
                    std::cout << "ac11 " << i + 1 << " " << get_cord_y(gap_str_cord) << " " << get_cord_y(pre_cord_start) << " " << endPosition(rstr[ri_pre + 1]) << "\n";
                    insertNewBamRecord(bam_records, g_id, bam_start, bam_strand); 
                }      
            }
            else if (flag & 1)
            {
               //NONE 
            }
            else if (flag & 2)
            {
                clip_segs(rstr[ri_pre], rstr[ri_pre + 1], 
                          pre_cord_start, _gap_parm, 1); 
                gap_end_cord = shift_cord(pre_cord_start, 
                                          endPosition(rstr[ri_pre]), 
                                          endPosition(rstr[ri_pre + 1]));
                std::cout << "ac3 " << get_cord_y(gap_end_cord) << "\n";
                if(insertGaps(gaps, gap_str_cord, gap_end_cord,
                              //rstr[ri_gap_h], 
                              //rstr[ri_gap_h + 1],
                              //rstr[ri_pre],
                              //rstr[ri_pre + 1],
                              length(bam_records) - 1,
                              thd_merge_gap))
                {
                    insertNewBamRecord(bam_records, g_id, bam_start, bam_strand); 
                }               
            }
        }
        else if (flag_pre & 2)
        {
            if (!flag)
            {
                clip_segs(rstr[ri_pre], rstr[ri_pre + 1], 
                          pre_cord_start, _gap_parm, -1);
                bam_start = get_cord_x(pre_cord_start) 
                         + beginPosition(rstr[ri_pre]);
                insertNewBamRecord(bam_records, 
                                   rstr[ri_pre], 
                                   rstr[ri_pre + 1],
                                   g_id, bam_start, bam_strand);

            }
            else if (flag & 1)
            {
                clip_segs(rstr[ri_pre], rstr[ri_pre + 1], 
                          pre_cord_start, _gap_parm, -1);     
                //detach_copy_row(rstr[ri_gap_h], rstr[ri_pre]);
                //detach_copy_row(rstr[ri_gap_h + 1], rstr[ri_pre + 1]); 
                gap_str_cord = shift_cord(pre_cord_start,
                                          beginPosition(rstr[ri_pre]),
                                          beginPosition(rstr[ri_pre + 1]));
                bam_start = get_cord_x(pre_cord_start) 
                         + beginPosition(rstr[ri_pre]);
                insertNewBamRecord(bam_records, 
                                   g_id, bam_start, bam_strand);
            }
            else if (flag & 2)
            {
                clip_segs(rstr[ri_pre], rstr[ri_pre + 1], 
                          pre_cord_start, _gap_parm, 0); 
                bam_start = get_cord_x(pre_cord_start) 
                         + beginPosition(rstr[ri_pre]);
                insertNewBamRecord(bam_records, 
                                   rstr[ri_pre], 
                                   rstr[ri_pre + 1],
                                   g_id, bam_start, bam_strand);  
            }
        }
        pre_cord_start = cord_start;
        pre_cord_end = cord_end;
        std::swap (ri, ri_pre); //swap the current and pre row id in the aligner.

    }
    printGaps(gaps.c_pairs);
    Score<int> score_scheme;
    int thd_clip_score = 80;
    int thd_reject_score = 130;
    int thd_accept_score = 140;
    int thd_accept_density = 16;
    align_gaps(bam_records, gaps, genomes, read, comrevRead, score_scheme, _gap_parm); 
    printCigarSrcLen(bam_records, "pscr_gaps1 ");
    return 0;
}

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
 uint64_t clip_window (String<Dna5> & genome,
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
    uint64_t returnCord = create_cord(genomeId, genomeStart + clip_ref, readStart + clip_read, strand);
    return returnCord;
}
