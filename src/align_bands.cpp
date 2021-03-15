#include <limits>
#include "base.h"
#include "cords.h"
#include "align_bands.h"

MergeCordsBandsParm::MergeCordsBandsParm()
{
    thd_search_depth = 10;
}

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
    uint64_t l_upper = l - band_upper;
    uint64_t l_lower = l - band_lower;
    return l * l - l_upper * l_upper / 2 - l_lower * l_lower / 2;
}
/*
 * Any band
 */
uint64_t calBandSize(uint64_t cord_str, uint64_t cord_end,
                     uint64_t band_lower, uint64_t band_upper)
{
    uint64_t l_x = get_cord_x(cord_end) - get_cord_x(cord_str);
    uint64_t l_y = get_cord_y(cord_end) - get_cord_y(cord_str);
    uint64_t l_x_upper = l_x - band_upper;
    uint64_t l_y_upper = l_y - band_upper;
    uint64_t l_x_lower = l_x - band_lower;
    uint64_t l_y_lower = l_y - band_lower;
    return l_x * l_y - l_x_upper * l_x_lower / 2 - l_y_upper * l_y_lower / 2;
}
/*
 * Sum of bands to be aligned
 */
uint64_t calBandsSize(String<uint64_t> & cords_str, String<uint64_t> & cords_end,
                      String<uint64_t> & bands_lower, String<uint64_t> & bands_upper)
{
    int size = 0;
    for (int i = 0; i < length(cords_str); i++)
    {
        size += calBandSize(cords_str[i], cords_end[i], bands_lower[i], bands_upper[i]);
    }
    return size;
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
/*
struct BandRecord
{

};

struct BandRecords
{
    String<_BandRecord> records; 

};
*/
/*
 * Only called by createAlignCords.
   Merge the cords for alignment if anchors are close or cords are close enough
   S1: [@cord11, @cord12):@band1 merge with
   S2: [@cord21, @cord22):@band2
   S1, S2 are required to be squares.
 * @thd_band_bound is the bound of the merged cords
   lower_bound upper_bound are supposed to be equal
 */
/*
int mergeCordsBands2(String<uint64_t> & cords_str,
                    String<uint64_t> & cords_end,
                    String<uint64_t> & bands_lower,
                    String<uint64_t> & bands_upper,
                    unsigned str_i, unsigned end_i,
                    unsigned thd_search_depth)
{
    BandRecords band_records(thd_search_depth);
    int min_area = std::numeric_limits<int>::max();
    for (unsigned i = str_i + 1; i < end_i; i++)
    {
        if (!isDiffCordsStrand(cords_str[i - 1], cords_str[i]) && 
            !_DefaultCord.isBlockEnd(cords_str[i - 1]))
        {
            BandRecord tmp_band_record(cords_str[i], cords_end[i], bands_lower[i], bands_upper[i]);
            BandRecord new_band_record(cords_str[i], cords_end[i], bands_lower[i], bands_upper[i]);
            int len = length(band_records.records);
            for (unsigned j = len - 1; j >= 0; j--)
            {
                int ii = i - len + j;
                BandRecord current_band_record(cords_str[ii], cords_end[ii], 
                    bands_lower[ii], bands_upper[ii]);
                BandRecord tmp_band_record;
                if (mergeCordBand(tmp_band_record, tmp_band_record, tmp_band_record))
                int area = tmp_band_record.band_area + band_records.records[j].band_area;    
                if (area < mini_area)
                {

                }
            }
        }
        else
        {
            band_records.clear();
            band_records.init();
        }
    }
    return 0 
}
*/
int mergeCordsBands1(String<uint64_t> & cords_str,
                     String<uint64_t> & cords_end,
                     String<uint64_t> & bands_lower,
                     String<uint64_t> & bands_upper,
                     unsigned str_i, unsigned end_i)
{
    return 0;
}
int mergeCordsBands(String<uint64_t> & cords_str,
                    String<uint64_t> & cords_end,
                    String<uint64_t> & bands_lower,
                    String<uint64_t> & bands_upper,
                    MergeCordsBandsParm & parm)
{
    int i_str = 1;
    for (int i = 2; i < length(cords_str); i++)
    {
        if (_DefaultCord.isBlockEnd(cords_str[i]) || isDiffCordsStrand(cords_str[i], cords_str[i - 1]))
        {
            mergeCordsBands1(cords_str, cords_end, bands_lower, bands_upper, i_str, i);
        }
        i_str = i + 1;
    }
    return 0;
}