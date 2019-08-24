#include <seqan/arg_parse.h>
#include <iostream>
#include <fstream>
#include <ctime>
#include "cords.h"
#include "pmpfinder.h"
#include "gap.h"
#include "align_interface.h"
#include "mapper.h"
#include "test_units.h"
//
//Check bamrecord cigar on base level
int check_cigar(StringSet<String<Dna5> > & genomes,
                 String<Dna5> & read, 
                 String<Dna5> & comrevRead,
                 String<uint64_t> & cords, //raw cords
                 String<BamAlignmentRecordLink> & bam_records)
{
    uint64_t str1, str2;
    String<Dna5> infix1;
    String<Dna5> infix2;
    int count_mat = 0;
    int count_mis = 0;
    int len = length(read) - 1;
    int f_new = 1;
    int j = 0;
    int it1 = 0;
    int it2 = 0;
    int base = 1;
    int f_b_end = 1;
    int seg_str_r;
    typedef std::pair<int, int> PairType ;
    String<PairType> covs;
    appendValue (covs, PairType(0,0));
    for (int i = 0; i < length(bam_records); i++)
    {
        std::cout << "chxxxb1 " << i<< " " << bam_records[i].isEnd() << " ";
        if (!(length(bam_records[i].cigar) == 0 ||
          ((length(bam_records[i].cigar) == 1) && 
                (bam_records[i].cigar[0].operation == 'S' || 
                bam_records[i].cigar[0].operation == 'H'))))
        {
            if (f_new)
            {
                count_mat = count_mis = 0;
                str1 = bam_records[i].beginPos;
                infix1 = infix(genomes[0], str1, length(genomes[0]));
                std::cout << "chxxb1" << str1<< "\n";
                if (bam_records[i].cigar[0].operation == 'S') 
                {
                    str2 = bam_records[i].cigar[0].count;
                    j = 1;
                }
                else
                {
                    str2 = 0;
                    j = 0;
                }
                if (bam_records[i].flag & 16) // - strand
                {
                    infix2 = infix(comrevRead, str2, length(read) - 1);
                }
                else
                {
                    infix2 = infix(read, str2, length(read) - 1);
                } 
                //std::cout << "nm0 " << infix1 << "\n";
                //std::cout << "nm1 " << infix2 << "\n";
                it1 = 0;
                it2 = 0;
                base = 1;
            }
            else
            {
                j = 0;
            }
            seg_str_r = it2;
            for (; j < length(bam_records[i].cigar); j++)
            {
                char op = bam_records[i].cigar[j].operation;
                int cnt = bam_records[i].cigar[j].count;
                if ((it1 + cnt >= length(infix1) && (op != 'I'))|| 
                    (it2 + cnt >= length(infix2) && (op != 'D')))
                {
                    dout << "bk1" << i << j << it1 << it2 << cnt << length(infix2) << "\n";
                    break;
                }
                if (op == 'D')
                {
                    for (int k = 0; k < cnt; k++)
                    std::cout << "chxx " << it1 + k << " " << infix1[it1 + k] << " | " << cnt << " " << op << "\n"; 
                    it1 += cnt;
                }
                else if (op == 'I')
                {
                    for (int k = 0; k < cnt; k++)
                    std::cout << "chxx " << it1 + k << " | " << infix2[it2 + k] << " " << cnt << " " << op  << "\n"; 
                    it2 += cnt;
                }
                else if (op == 'X')
                {
                    for (int k = 0; k < cnt; k++)
                    std::cout << "chxx " << it1 + k << " " << infix1[it1 + k] << " " << infix2[it2 + k] << " " << cnt << " " << op  << "\n"; 
                    it1 += cnt;
                    it2 += cnt;
                }
                else if (op == '=')
                {
                    for (int k = 0; k < cnt; k++)
                    {
                        std::cout << "chxx " << it1 << " " << infix1[it1] << " " << infix2[it2] << " " << cnt << " " << op << " "; 
                        if (infix1[it1] != infix2[it2]) 
                        {
                            ++count_mis;
                            std::cout << "xxxxxxxxxxxxxxxxxxxxxx\n";
                        }
                        else
                        {
                            ++count_mat;
                            std::cout << "\n";
                        }
                        it1++;
                        it2++;
                    }
                }
            }
        }
        if (bam_records[i].isEnd())
        {
            if (count_mis != 0)
            {
                dout << "cb1 " << i << " " << str1 << str2 << " " << it1 << " " << it2 << length(infix1) << length(infix2) << count_mat << " " << count_mis << " " << "\n";
                std::cerr << "cb1\n";
            }
            else
            {
                dout << "cb3 " << i << " " << str1 << " " << str2 << it1 << " " << it2 << length(infix1) << length(infix2) << count_mat << " " << count_mis << " " << "\n";
                appendValue (covs, PairType(seg_str_r, it2));
            }
            f_new = 1;
        }
        else 
        {
            dout << "chxxb0 " << i << " " << str1 << " " << str2 << " " << bam_records[i + 1].beginPos << seg_str_r << it2 << "\n";
     //       appendValue (covs, PairType(seg_str_r, it2));
            i = bam_records[i].next() - 1;
            f_new = 0;
        }
    }
    //appendValue (covs, PairType (length(read), length(read)));
    //std::sort (begin(covs), end(covs), [](PairType a, PairType b){return a.first < b.first;});
    int covs_s = 0;
    /*
    for (int i = 0; i < length(covs) - 1; i++)
    {
        int d = covs[i + 1].first - covs[i].second;
        if (d > 0)
        {
            covs_s += d;
        }
    }
    */
    dout << "c_s" << float(covs_s) / length(read) << "\n";
}
/*
int test_clip_anchors()
{

}
*/
int check_index1(StringSet<String<Dna5> > & seqs, IndexDynamic & index, unsigned threads)
{
    HIndex hindex = index.hindex;
    LShape shape(hindex.getShape());
    dout << "lensq" << length(seqs[0]) << "\n";
    uint thd_step = 1;
    createHIndex(seqs, hindex, 0, length(seqs), threads, thd_step, false);

    for (auto seq : seqs)
    {
        int64_t count = 0;
        hashInit (shape, begin(seq));
        for (int i = 0; i < length(seq) - shape.span; i++)
        {
            hashNexth (shape, begin(seq) + i);
            if (++count > thd_step)
            {
                hashNextX (shape, begin(seq) + i);
                uint f_find = 0;
                uint64_t pos = getXDir(hindex, shape.XValue, shape.YValue);
                uint64_t ptr = _DefaultHs.getHeadPtr(hindex.ysa[pos-1]);
                while ((_DefaultHs.getHsBodyY(hindex.ysa[pos]) == shape.YValue || 
                        _DefaultHs.getHsBodyY(hindex.ysa[pos]) == 0))
                {
                    if (((hindex.ysa[pos] & _DefaultHsBase.bodyCodeFlag) >>_DefaultHsBase.bodyCodeBit) ^ shape.strand)
                    {

                    }
                    else
                    {    
                        f_find = 1;
                        dout << "f_find" << f_find << "\n";
                    }
                    if (++pos > length(hindex.ysa) - 1)
                    {
                        break;
                    }
                }
                if (f_find == 0)
                {
                    dout << "cidx2_1" << i  << "\n";
                    //return 1;
                }
                count = 0;
            }
        }
    }
    return 0;

}

/* test Dynamic index2 (direct index)

*/
int check_index2(StringSet<String<Dna5> > & seqs, IndexDynamic & index, unsigned threads)
{
    int64_t thd_min_step = 0;
    int64_t thd_max_step = 0;
    int64_t thd_omit_block = 1ULL << 62; //max infi
    DIndex dindex = index.dindex;
    LShape shape(dindex.getShape());
    dout << "lensq" << length(seqs[0]) << "\n";

    createDIndex(seqs, dindex, thd_min_step, thd_max_step, thd_omit_block, 0, length(seqs), threads);

    dout << "lenseq1" << back(dindex.getDir()) << "\n";

    for (auto seq : seqs)
    {
        int64_t count = 0;
        hashInit (shape, begin(seq));
        for (int i = 0; i < length(seq) - shape.span; i++)
        {
            hashNexth (shape, begin(seq) + i);
            if (++count > thd_min_step)
            {
                hashNextX (shape, begin(seq) + i);
                //dout << "shape" << shape.XValue << "\n";
                int64_t str_ = queryHsStr(dindex, shape.XValue);
                int64_t end_ = queryHsEnd(dindex, shape.XValue);
                int f_find = 0;
                //dout << "x4xxxxxxxxxx " << i << end_ -  str_ << "\n";
                for (int j = str_; j < end_; j++)
                {
                    int x = get_cord_x (dindex.getHs()[j]);
                    //dout << "x4" << x << j << "\n";
                    if (x == i)
                    {
                        f_find = 1;
                        //dout << "find " << i << x << "\n";
                    }
                }
                if (f_find == 0)
                {
                    dout << "cidx2_1" << i << end_ - str_ << "\n";
                    return 1;
                }
                count = 0;
            }
        }
    }
    return 0;

}
