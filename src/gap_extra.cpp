#include "gap_util.h"
/********************** START: Extend mapping for DUP *************************/
/*
 * @gap_str, @gap_end are required to have the same strand
 */
uint64_t try_dup_filter_(String<uint64_t> & dup_tiles, uint64_t gap_str, uint64_t gap_end, 
    int direction)
{
    int64_t thd_anchor_err = 80; //todo;change this
    int64_t anchor = direction < 0 ? get_tile_x(gap_str) - get_tile_y(gap_str) : get_tile_x(gap_end) - get_tile_y(gap_end);
    uint64_t ii = 0;
    for (uint64_t i = 0; i < length(dup_tiles); i++)
    {
        int64_t dup_anchor = get_tile_x(dup_tiles[i]) - get_tile_y(dup_tiles[i]);
        if ((!get_tile_strand(dup_tiles[i] ^ gap_str)) && (std::abs(dup_anchor- anchor) > thd_anchor_err))
        {
            if (is_tile_end(dup_tiles[i]) && i > ii)
            {
                set_tile_end(dup_tiles[i - ii - 1]);
            }
            ++ii;
        }
        else
        {
            dup_tiles[i - ii] = dup_tiles[i];
        }
    }
    resize (dup_tiles, length(dup_tiles) - ii);
    return 0;
}
/**
 * Additional duplication (only) mapping in [gap_str,) towards right and (,gap_end] towards right.
 */
int try_dup(String<Dna5> & seq,
            String<Dna5> & read,
            String<Dna5> & comstr,
            StringSet<FeaturesDynamic > & f1,
            StringSet<FeaturesDynamic > & f2,
            String<uint64_t> & tiles_left, //result
            String<uint64_t> & tiles_rght, //result
            uint64_t gap_str,
            uint64_t gap_end,
            float thd_err_rate,
            uint64_t thd_tile_size,
            GapParms & gap_parms)
            //float band_ratio,
            //int direction)
{
    int64_t dy;
    if (get_tile_strand(gap_str ^ gap_end))
    {
        dy = length(read) - 1 - get_tile_y(gap_end) - get_tile_y(gap_str);
    }
    else
    {
        dy = get_tile_y(gap_end) - get_tile_y(gap_str);
    }
    if (dy < 0)
    {
        return 1;
    }
    int64_t shift_y = dy; 
    int64_t shift_x = dy * (1 + thd_err_rate);

    //extend tile1 towards right
    uint64_t dup_str = gap_str;
    uint64_t dup_end = shift_tile (gap_str, shift_x, shift_y);
    int direction = 1;
    mapInterval(seq, read, comstr,  tiles_left, f1, f2, 
                 dup_str, dup_end, LLMIN, LLMAX, direction,
                 gap_parms);
    if (!empty(tiles_left))
    {
        set_tile_end(back(tiles_left));
    }
    try_dup_filter_(tiles_left, gap_str, gap_end, -1);

    //extend tile2 towards left
    dup_end = gap_end;
    dup_str = shift_tile(gap_end, -shift_x, -shift_y);
    direction = -1;
    String<uint64_t> tiles2;
    mapInterval(seq, read, comstr, tiles_rght, f1, f2,
                 dup_str, dup_end,LLMIN, LLMAX, direction,
                 gap_parms);
    if (!empty(tiles_rght))
    {
        remove_tile_sgn_end(back(tiles_rght));
    }
    try_dup_filter_(tiles_rght, gap_str, gap_end, 1);
    (void)thd_tile_size;
    return 0;
}

/********************** START: Extend mapping for INV *************************/
/*
int __extendsIntervalClipOverlapsInv(
        String<uint64_t> & chain1, String<uint64_t> & chain2, String<uint64_t> & chain3,
        String<uint64_t> & chain4, int shape_len, bool f_clip, uint64_t(*getX)(uint64_t),
        uint64_t(*getY)(uint64_t), uint64_t read_len, GapParms & gap_parms)
{
    String<int> gaps_score11, gaps_score12;
    String<int> gaps_score21, gaps_score22;
    String<int> gaps_score31, gaps_score32;
    String<int> gaps_score41, gaps_score42;
    accumulateSimpleGapScore1(chain1, gaps_score11, shape_len, getX, gap_parms);
    accumulateSimpleGapScore1(chain1, gaps_score12, shape_len, getY, gap_parms);
    accumulateSimpleGapScore1(chain2, gaps_score21, shape_len, getX, gap_parms);
    accumulateSimpleGapScore1(chain2, gaps_score22, shape_len, getY, gap_parms);
    accumulateSimpleGapScore1(chain3, gaps_score31, shape_len, getX, gap_parms);
    accumulateSimpleGapScore1(chain3, gaps_score32, shape_len, getY, gap_parms);
    accumulateSimpleGapScore1(chain4, gaps_score41, shape_len, getX, gap_parms);
    accumulateSimpleGapScore1(chain4, gaps_score42, shape_len, getY, gap_parms);
    int i1 = i2 = i3 = i4 = 0;
    int max_score = 0;
    int is[4]={-1,-1,-1,-1};
    for (; i1 < length(chain1); i1++) 
    {
        uint64_t x1 = getX(chain1[i1]);
        uint64_t y1 = getY(chain1[i1]);

        uint64_t x2_lower = x1;   
        uint64_t x2_upper = x2_lower + gap_parms.thd_eicos_clip_dxy;
        uint64_t y3_upper = read_len - 1 - y1; 
        uint64_t y3_lower = y3_upper - gap_parms.thd_eicos_clip_dxy;

        score1 = gaps_score11[i1] + gaps_score12[i1];
        for (i2; i2 < length(chain2); j++)
        {
            x2 = getX(chain2[j]);
            y2 = getY(chain2[j]);
            if (x2 < x2_lower)
            {
                continue;
            }
            else if (x2 > x2_upper)
            {
                break;
            }
            uint64_t y4_lower = read_len - y2 - 1;
            uint64_t y4_upper = y4_lower + gap_parms.thd_eicos_clip_dxy;
            score2 = back(gaps_score21) - gaps_score21[i2] + 
                     back(gaps_score22) - gaps_score22[i2];

            for (i4; i4 < length(chain4); j++)
            {
                x4 = getX(chain4[i4]);
                y4 = getY(chain4[i4]);
                if (y4 < y4_lower)
                {
                    continue;
                }
                else if (y4 > y4_upper)
                {
                    break;
                }
                uint64_t x3_upper = x4;
                uint64_t x3_lower = x3_upper - gap_parms.thd_eicos_clip_dxy;
                score4 = back(gaps_score41) - gaps_score41[i4] +
                         back(gaps_score42) - gaps_score42[i4];
                for (i3; i3 < length(chain3); i3++)
                {
                    if (x3 < x3_lower || y3 < y3_lower)
                    {
                        continue;
                    }
                    if (x3 > x3_upper || y3 > y3_upper)
                    {
                        break;
                    }
                    score3 = gaps_score31[i3] + gaps_score32[i3];
                    score = score1 + score2 + score3 + score4;
                    if(score > max_score) 
                    {
                        max_score = score;
                        is[0] = i1;
                        is[1] = i2;
                        is[2] = i3;
                        is[4] = i4;
                    }
                }
            }
        }
    }
    if (i[0] >= 0)
    {
        if (f_clip)
        {
            erase(chain1, is[0] + 1, length(chains1));
            erase(chain2, 0, is[1] + 1);
            erase(chain3, is[2] + 1, length(chains3));
            erase(chain4, 0, is[3] + 1);
        }
        return 0;
    } 
    else
    {
        return -1;
    }
}
int extendsIntervalClipOverlapsInv_(
        String<uint64_t> & chain1, String<uint64_t> & chain2, String<uint64_t> & chain3,
        String<uint64_t> & chain4, int shape_len, bool f_clip, uint64_t(*getX)(uint64_t),
        uint64_t(*getY)(uint64_t), GapParms & gap_parms)
{
 
    clipChain (chain1, shape_len, g_map_left, true, getX, getY, gap_parms);
    clipChain (chain2, shape_len, g_map_left, true, getX, getY, gap_parms);
    clipChain (chain3, shape_len, g_map_left, true, getX, getY, gap_parms);
    clipChain (chain4, shape_len, g_map_left, true, getX, getY, gap_parms);
    return 0;
}
int extendsIntervalMapOverlapsInv_(
        String<uint64_t> & chain1, String<uint64_t> & chain2, String<uint64_t> & chain3,
        String<uint64_t> & chain4, int shape_len, bool f_clip, uint64_t(*getX)(uint64_t),
        uint64_t(*getY)(uint64_t), GapParms & gap_parms)
{

}
*/
/*
 * A simple class to record start and end of invs
 *
class InvStack
{
    String<int> inv_stack; 
public:
    int match(String<uint64_t> & cords, int i, uint64_t strand1, uint64_t strand2)
    {
        if (!(strand1 ^ strnad2))
        {
            return 0;
        }
        if (empty(inv_stack))
        {
            strand1 = strand1 ? uint(1) : uint(0);
            strand2 = strand1 ? uint(1) : uint(0);
            appendValue(inv_stack, strand1 + (strand2 << 1) + (i << 2));
            return 0;
        }
        else 
        {
            uint64_t record = back(inv_stack);
            str_strand1 = record & 1; 
            str_strand2 = (record >> 1)& 1;
            str_i = (record >> 2); 
            if ((str_strand1 ^ strand1) && (str_strand2 ^ strand2) )
            {
                return str_i;
            }
            else
            {
                appendValue(inv_stack, strand1 + (strand2 << 1) + (i << 2));
                return 0;
            }
        }
        return 0;
    };
    void clear()
    {
        clear(inv_stack);
    };
    void update(int i, int delta) 
    {
        for (int j = 0; j < inv_stack; j++)
        {
            if ((inv_stack[j] >> 2) >= i)
            {
                inv_stack[j] += (delta << 2);
            }
        }
    }
}
int mapGapGlobal(String<uint64_t> & tiles1,
                 String<uint64_t> & tiles2,
                 String<uint64_t> & tiles3,
                 String<uint64_t> & tiles4,
                 uint64_t cord_str1, uint64_t cord_end1, 
                 uint64_t cord_str2, uint64_t cord_end2,
                 uint64_t cord_str3, uint64_t cord_end3, 
                 uint64_t cord_str4, uint64_t cord_end4,
                 uint64_t read_len,
                 GapParms & gap_parms
                 )
{
    InvStack inv_stack;
    for (unsigned i = 1; i < length(cords); i++)     
    {
        if (_DefaultCord.isBlockEnd(cords[i - 1]))
        {
            inv_stack.clear();
        }
        uint64_t strand1 = get_cord_strand(cords[i - 1]);
        uint64_t strand2 = get_cord_strand(cords[i]);
        if (get_cord_strand(cords[i - 1] ^ cords[i]))
        {
            int inv_str = inv_stack.match(i, strand1, strand2);
            if (inv_str)
            {
                extendsIntervalClipOverlapsInv_();
                if (!empty())
                {
                }
            }
        }
    }
    return 0;
}
*/