#include "base.h"
#include "cords.h"
#include "align_bands.h"

/*--------------------  Merge bands of cords for alignment  --------------------
  For simplicity and efficency, only merge bands of 45 degree.
 */
/*
 * Caculate the size of bands of 45 degree
 */
uint64_t calBandSize45(uint64_t cord_str, uint64_t cord_end,
                       uint64_t band_lower, uint64_t  band_upper)
{
    uint64_t l = get_cord_x(cord_end) - get_cord_x(cord_str);
    return l * l - (l - band_upper) * (l - band_upper) / 2 
                 - (l - band_lower) * (l - band_lower) / 2;
}

/*
 * Merge region S1=[@cord_str1, @cord_end1) and S2=[@cord_str2, @cord_end2)
 * @band11, @band12 are lower and upper bands of S1 
 * @band21, @band22 are lower and upper bands of S2 
 * The diagnal of merged cord is required to be < d
 */
int mergeCordBand (uint64_t & cord_str1, uint64_t & cord_end1,
                   uint64_t & cord_str2, uint64_t & cord_end2,
                   int & band_lower1, int & band_upper1, 
                   int & band_lower2, int & band_upper2, 
                   int thd_max_band)
{
    return 0;
}


int mergeCordsBands(String<uint64_t> & cords_str,
               String<uint64_t> & cords_end,
               String<int> & bands1,
               String<int> & bands2,
               int thd_max_band)
{
    resize(bands1, length(cords_str))
    resize(bands2, length(cords_end))
    int ii = 0;
    for (int i = 0; i < length(cords_str); i++)
    {
        if(mergeCordBand(cords_str[ii], cords_end[ii], cords_str[i], cords_ends[i], 
                bands1[ii], bands2[ii], bands1[i], bands2[i], thd_max_band))
        {
            ++ii;
        }
    }
    resize(cords_str, ii);
    resize(cords_end, ii);
    resize(bands1, ii);
    resize(bands2, ii);
    return 0;
}


class __MergeBandsBuffer
{
    String<std::pair<unsigned, unsigned> > buffer; 
    initBuffer();
}
/*
 * Only called by createAlignCords.
   Merge the cords for alignment if anchors are close or cords are close enough
   S1: [@cord11, @cord12):@band1 merge with
   S2: [@cord21, @cord22):@band2
   S1, S2 are required to be squares.
 * @thd_band_bound is the bound of the merged cords
   lower_bound upper_bound are supposed to be equal
 */
int mergeAlignCords(String<uint64_t> & cords_str1,
                    String<uint64_t> & cords_end1,
                    String<uint64_t> & cords_str2,
                    String<uint64_t> & cords_end2,
                    String<uint64_t> & bands_lower1,
                    String<uint64_t> & bands_upper1,
                    String<uint64_t> & bands_lower2,
                    String<uint64_t> & bands_upper2,
                    unsigned str_i, unsigned end_i,
                    unsigned thd_search_depth)
{
    unsigned block_str = str_i;
    String<std::pair<unsigned, unsigned> > buffer;
    resize(buffer, thd_search_depth);
    for (unsigned i = str_i + 1; i < end_i; i++)
    {
        if (!isDiffCordsStrand(cords_str1[i - 1], cords_str1[i]) && !isBlockEnd(cords_str1[i - 1]))
        {
            for (unsigned j = buffer.len() - 1; j > 0; j--)
            {
                if (createMergedBands(cords_str[j - 1], cords_end[j - 1],
                                      cords_str[j]), cords_end[j], 
                                      band_lower, band_upper)
                {
                    
                }
            }
        }
        else
        {
            //clearBestBuffer()
            //initBestBuffer()
        }
    }
    return 0 
}
