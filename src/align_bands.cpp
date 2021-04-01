#include <limits>
#include "base.h"
#include "cords.h"
#include "align_bands.h"

struct LineSegment
{
    uint64_t str_x;
    uint64_t str_y;
    uint64_t end_x;
    uint64_t end_y;
    LineSegment(uint64_t x1, uint64_t y1, uint64_t x2, uint64_t y2)
    {
        str_x = x1;
        str_y = y1;
        end_x = x2;
        end_y = y2;
    }
};

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
    dout << "cbz1" << l_x << l_y << l_x_lower << l_y_lower << l_x_upper << l_y_upper << "\n";
    return l_x * l_y - (l_x_lower * l_y_lower  + l_x_upper * l_y_upper) / 2;
}
/*
 * Sum of bands to be aligned
 */
uint64_t calBandsSize(String<uint64_t> & cords_str, String<uint64_t> & cords_end,
                      String<uint64_t> & bands_lower, String<uint64_t> & bands_upper,
                      unsigned i_str, unsigned i_end)
{
    uint64_t size = 0;
    for (unsigned i = i_str; i < i_end; i++)
    {
        uint64_t new_size = calBandSize(cords_str[i], cords_end[i], bands_lower[i], bands_upper[i]);
        size += new_size;
        dout << "cbsz1" << new_size << "\n";
    }
    return size;
}
int isColinear(uint64_t x1, uint64_t y1,
               uint64_t x2, uint64_t y2,
               uint64_t x3, uint64_t y3)
{
    return (x1 - x2) * (y2 - y3) == (x2 - x3) * (y1 - y2);
}
int isColinear(uint64_t x11, uint64_t y11,
               uint64_t x12, uint64_t y12,
               uint64_t x21, uint64_t y21,
               uint64_t x22, uint64_t y22)
{
    return isColinear(x11, y11, x12, y12, x21, y21) &&
           isColinear(x12, y12, x21, y21, x22, y22);
}
int isLineSegmentColinear(LineSegment & line1, LineSegment & line2)
{
    return isColinear(line1.str_x, line1.str_y, line1.end_x, line1.end_y,
                      line2.str_x, line2.str_y, line2.end_x, line2.end_y);
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
/*
 * Shortcut function:
 * Given the start and end of the rectangle, get the start and end of the band within
   the rectangle
 */
LineSegment getRectangleBandStrEnd(uint64_t x_str, uint64_t y_str,
        uint64_t x_end, uint64_t y_end, uint64_t band, uint64_t f_l_u)
{
    if (f_l_u < 0) //lower band
    {
        return  LineSegment(x_str, y_str + band, x_end - band, y_end);
    }
    else
    {
        return  LineSegment(x_str + band, y_str, x_end, y_end - band);
    }
}
int mergeCordsBands1(String<uint64_t> & cords_str,
                     String<uint64_t> & cords_end,
                     String<uint64_t> & bands_lower,
                     String<uint64_t> & bands_upper,
                     unsigned i_str, unsigned i_end)
{
    unsigned it = i_str;
    uint64_t x11 = get_cord_x(cords_str[i_str]); //predecessor after merge
    uint64_t y11 = get_cord_y(cords_str[i_str]);
    uint64_t x12 = get_cord_x(cords_end[i_str]);
    uint64_t y12 = get_cord_y(cords_end[i_str]);
    uint64_t x21 = 0;                          //current
    uint64_t y21 = 0;
    uint64_t x22 = 0;
    uint64_t y22 = 0;

    //<<debug
    int area1 = calBandsSize (cords_str, cords_end, bands_lower, bands_upper, i_str, i_end);
    //>>debug
    String<unsigned> del_list; //temp array to record ith item to delete
    for (unsigned i = i_str + 1; i < i_end; i++)
    {
        x21 = get_cord_x(cords_str[i]);
        y21 = get_cord_y(cords_str[i]);
        x22 = get_cord_x(cords_end[i]);
        y22 = get_cord_y(cords_end[i]);
        dout << "mcb10" << x12 << x21 << y12 << y21 << "\n";

        LineSegment band_lower1 = 
            getRectangleBandStrEnd(x11, y11, x12, y12, bands_lower[it], -1);
        LineSegment band_upper1 = 
            getRectangleBandStrEnd(x11, y11, x12, y12, bands_upper[it],  1);
        LineSegment band_lower2 = 
            getRectangleBandStrEnd(x21, y21, x22, y22, bands_lower[i], -1);
        LineSegment band_upper2 = 
            getRectangleBandStrEnd(x21, y21, x22, y22, bands_upper[i],  1);
       
        if(x12 > x21 && y12 > y21 &&  //when current overlaps with last block, thus skip ins,del 
           isLineSegmentColinear(band_lower1, band_lower2) &&
           isLineSegmentColinear(band_upper1, band_upper2))
        {
            cords_end[it] = cords_end[i];
            //appendValue(del_list, i);
            dout << "mcoline" << x11 << y11 << x12 << y12 << x21 << y21 << x22 << y22 << float(y12-y11)/(x12-x11) << float(y21 - y12)/(x21 - x12) << float(y22 - y21)/(x22 - x21) << x12 - x21 << y12 - y21 << "\n";
            x12 = x22;
            y12 = y22;
        }
        else 
        {
            /*
            if (x12 - x21 > (x12 - x11) * thd_max_overlap_x &&
                y12 - y21 > (y12 - y11) * thd_max_overlap_y)
            {

            }
            */
            ++it;
            cords_str[it] = cords_str[i];
            cords_end[it] = cords_end[i];
            bands_lower[it] = bands_lower[i];
            bands_upper[it] = bands_upper[i];
            x11 = x21;
            y11 = y21;
            x12 = x22;
            y12 = y22;
        }
        dout << "mcb11" << x11 << y11 << x21 << y21 << x21 - x11 << y21 - y11 << bands_lower[i] << bands_upper[i] << "\n";
    }
    erase(cords_str, it + 1, i_end);
    erase(cords_end, it + 1, i_end);
    erase(bands_lower, it + 1, i_end);
    erase(bands_upper, it + 1, i_end);
    //<<debug
    print_cords(cords_str, "mcb13");
    int area2 = calBandsSize (cords_str, cords_end, bands_lower, bands_upper, i_str, i_str + it + 1);
    dout << "mcb12" << area1 << area2 << length(cords_str) << "\n";
    //>>debug
    return 0;
}
int mergeCordsBands(String<uint64_t> & cords_str,
                    String<uint64_t> & cords_end,
                    String<uint64_t> & bands_lower,
                    String<uint64_t> & bands_upper,
                    MergeCordsBandsParm & parm)
{
    int i_str = 1;
    for (unsigned i = 2; i < length(cords_str); i++)
    {
        if (i == length(cords_str) - 1 || _DefaultCord.isBlockEnd(cords_str[i]) || 
            (i < length(cords_str) - 1 && isDiffCordsStrand(cords_str[i], cords_str[i + 1])))
        {
            mergeCordsBands1(cords_str, cords_end, bands_lower, bands_upper, i_str, i + 1);
            dout << "mcb01" << i_str << i << "\n";
            i_str = i + 1;
        }
    }
    return 0;
}