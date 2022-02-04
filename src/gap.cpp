#include <utility> 
#include <seqan/align.h>
#include "base.h"
#include "shape_extend.h"
#include "cords.h"
#include "gap_util.h"
#include "gap_extra.h"
//#include "gap.h"

/*----------------------  Gap main  ----------------------*/
/*
 * Map and resolve the gap specified by [@gap_str, @gap_end)
 */
int mapGap_ (StringSet<String<Dna5> > & seqs, 
             String<Dna5> & read, 
             String<Dna5> & comstr,
             uint64_t gap_str, 
             uint64_t gap_end, 
             StringSet<FeaturesDynamic > & f1, 
             StringSet<FeaturesDynamic > & f2,
             String<uint64_t> & tiles_str, 
             String<uint64_t> & tiles_end, 
             String<uint64_t> & clips, //results 
             int direction,
             int64_t thd_dxy_min,
             GapParms & gap_parms)
{
    unused(clips);
    CmpInt64 g_cmpll;
    //float thd_da_zero = gap_parms.thd_err;
    //float thd_da_rate = 0.1;
    clear(tiles_str);
    clear(tiles_end);
    //dout << "mg3" << get_cord_y(gap_str) << get_cord_y(gap_end) <<direction << "\n";
    _DefaultHit.unsetBlockEnd(gap_str); //remove cord sgn, format cord to tiles
    _DefaultHit.unsetBlockEnd(gap_end);

    String <Dna5> & ref = seqs[get_cord_id(gap_str)];
    int64_t x1 = get_cord_x(gap_str);
    int64_t x2 = get_cord_x(gap_end);
    int64_t y1 = get_cord_y(gap_str);
    int64_t y2 = get_cord_y(gap_end);
    int64_t shift_x, shift_y;
    String<uint64_t> sp_tiles; 
    if (get_cord_strand(gap_str ^ gap_end))
    {
        if (direction != g_map_closed)
        {
            return -1; //this case is not allowed, 
            //since the gap_str or gap_end is infi cord
        }
        int64_t thd_max_extend1 = 500; //when x2 < x1
        int64_t thd_max_extend2 = 5000; //normal case 

        String<uint64_t> tiles_str1;
        String<uint64_t> tiles_str2;
        String<uint64_t> tiles_end1;
        String<uint64_t> tiles_end2;
        int direction1 = g_map_rght;
        int direction2 = g_map_left;

        shift_x = (x2 - x1 > 0) ? 
            std::min({thd_max_extend2, int64_t(length(ref) - 1 - get_cord_x(gap_str)), x2 - x1}) : thd_max_extend1; 
        g_cmpll.min(shift_y, (x2 - x1) * (1 + gap_parms.thd_err)) << int64_t(length(read) - 1 - get_cord_y(gap_str));
        shift_x = std::max(shift_x, int64_t(0));
        shift_y = std::max(shift_y, int64_t(0));
        uint64_t gap_str1 = gap_str;
        uint64_t gap_end1 = shift_cord (gap_str, shift_x, shift_y);
        //gap_parms.thd_gmsa_d_anchor_rate = 0.1;
        mapExtend (seqs, read, comstr, f1, f2, tiles_str1, tiles_end1, 
                    gap_str1, gap_end1, direction1, gap_parms);
        //g_cmpll.min(shift_x, x2 - x1) << int64_t(get_cord_x(gap_end)) << int64_t(5000);
        shift_x = (x2 - x1 > 0) ? std::min({x2 - x1, int64_t(get_cord_x(gap_end)), thd_max_extend2}) : thd_max_extend1;
        //g_cmpll.min(shift_x, x2 - x1) << int64_t(get_cord_x(gap_end));
        g_cmpll.min(shift_y, (x2 - x1) * (1 + gap_parms.thd_err)) << int64_t(get_cord_y(gap_end));
        shift_x = std::max(shift_x, int64_t(0));
        shift_y = std::max(shift_y, int64_t(0));
        uint64_t gap_str2 = shift_cord (gap_end, -shift_x, -shift_y);
        uint64_t gap_end2 = gap_end;
        mapExtend (seqs, read, comstr, f1, f2, tiles_str2, tiles_end2, 
                    gap_str2, gap_end2, direction2, gap_parms);

        if (!empty(tiles_str1))
        {
            append(tiles_str, tiles_str1);
            append(tiles_end, tiles_end1);
        }
        if (!empty(tiles_str2))
        {
            append(tiles_str, tiles_str2);
            append(tiles_end, tiles_end2);
        }
    }
    else if (y1 > (int64_t)length(read) - 1 || y2 > (int64_t)length(read) - 1 ||
             x1 > (int64_t)length(ref) - 1  || x2 > (int64_t)length(ref) - 1)
    {
        //return 1;
        //todo
    }
    else if (x1 < x2 && y1 < y2)
    {
        int64_t danc = x1 - x2 - y1 + y2;
        if (std::abs(danc) > gap_parms.thd_mg1_danc_indel &&
            direction == g_map_closed) //ins/del 
        {
            //!add declaration of variables here
            ChainScoreMetric chn_score1_tmp = gap_parms.chn_score1; //used to restore the original parms after mapExtend
            ChainScoreMetric chn_score2_tmp = gap_parms.chn_score2; //used to restore the original parms after mapExtend
            gap_parms.chn_score1.thd_min_chain_len = 1;
            gap_parms.chn_score1.thd_abort_score = 0;
            gap_parms.chn_score1.getScore = &getGapAnchorsChainScore2;
            gap_parms.chn_score2.thd_abort_score = 0;
            gap_parms.chn_score2.getScore2 = &getGapBlocksChainScore3;

            String<uint64_t> tiles_str1;
            String<uint64_t> tiles_str2;
            String<uint64_t> tiles_end1;
            String<uint64_t> tiles_end2;
            
            uint64_t gap_str1, gap_str2, gap_end1, gap_end2;
            if (danc > 0) //!ins
            {
                shift_y = std::min({std::max(y2 - y1, int64_t(0)),  gap_parms.thd_max_extend2, 
                    int64_t(length(read) - y1 - 1)});
                shift_x = std::min({int64_t(shift_y * (1 + gap_parms.thd_err)), 
                    gap_parms.thd_max_extend2, int64_t(length(ref) - x1 - 1)});
                gap_str1 = gap_str;
                gap_end1 = shift_cord (gap_str, shift_x, shift_y);
                shift_y = std::min({std::max(y2 - y1, int64_t(0)), gap_parms.thd_max_extend2, int64_t(y2)});
                shift_x = std::min({int64_t(shift_y * (1 + gap_parms.thd_err)), gap_parms.thd_max_extend2, int64_t(x2)});          
                gap_str2 = shift_cord (gap_end, -shift_x, -shift_y);
                gap_end2 = gap_end;
            }
            else //!del
            {
                shift_x = std::min({std::max(x2 - x1, int64_t(0)), gap_parms.thd_max_extend2, int64_t(length(ref) - x1 - 1)});
                shift_y = std::min({int64_t(shift_x * (1 + gap_parms.thd_err)), gap_parms.thd_max_extend2, int64_t(length(read) - y1 - 1)});
                gap_str1 = gap_str;
                gap_end1 = shift_cord (gap_str, shift_x, shift_y);
                shift_x = std::min({std::max(x2 - x1, int64_t(0)), gap_parms.thd_max_extend2, int64_t(x2)});
                shift_y = std::min({int64_t(shift_x * (1 + gap_parms.thd_err)), gap_parms.thd_max_extend2, int64_t(y2)});          
                gap_str2 = shift_cord (gap_end, -shift_x, -shift_y);
                gap_end2 = gap_end;
            }
            //!add map process here
            mapExtends (seqs, read, comstr, f1, f2,  
                        tiles_str1, tiles_end1, 
                        tiles_str2, tiles_end2, 
                        gap_str1, gap_end1,
                        gap_str2, gap_end2, 
                        thd_dxy_min, gap_parms);

            //!add post process here            
            if (!empty(tiles_str1))
            {
                append(tiles_str, tiles_str1);
                append(tiles_end, tiles_end1);
                remove_tile_sgn(back(tiles_str));
                remove_tile_sgn(back(tiles_end));

            }
            if (!empty(tiles_str2))
            {
                remove_tile_sgn(tiles_str2[0]);
                remove_tile_sgn(tiles_end2[0]);
                append(tiles_str, tiles_str2);
                append(tiles_end, tiles_end2);
            }
            gap_parms.chn_score2 = chn_score2_tmp;
            gap_parms.chn_score1 = chn_score1_tmp;
            //dout << "gap_indel" << get_cord_y(gap_str1) << "\n";
        }
        else 
        {
            //mapInterval(ref, read, comstr, tiles_str, f1, f2, gap_str, gap_end, LLMIN, LLMAX, direction, gap_parms);
            /*
            reform_tiles(ref, read, comstr, 
                         tiles_str, tiles_end,
                         clips, sp_tiles,
                         gap_str, gap_end, direction, 
                         thd_cord_gap, thd_tile_size, thd_err_rate);
                         */
            //try dup for each sp_tiles element.
            for (uint64_t i = 0; int(i) < getClipsLen(sp_tiles); i++)
            {
                
              //  String <uint64_t> dup_tiles_left_str;
              //  String <uint64_t> dup_tiles_rght_str;
              //  String <uint64_t> dup_tiles_left_end;
              //  String <uint64_t> dup_tiles_rght_end;
              //  String <uint64_t> dup_clips;  //not used
              //  String <uint64_t> dup_sp_tiles; //not used
              //  int dup_direction = g_map_closed;
              //  uint64_t it_dup_str = getClipStr(sp_tiles, i) + new_dups_count;
              //  uint64_t dup_str = tiles_str[getClipStr(sp_tiles, i) + new_dups_count];
              //  uint64_t dup_end = tiles_str[getClipEnd(sp_tiles, i) + new_dups_count];
              //  //if ()
              //  dup_end = shift_tile(dup_end, thd_tile_size, thd_tile_size);
              //  dout << "dpstrend" << getClipStr(sp_tiles, i) << getClipEnd(sp_tiles, i) << new_dups_count << "\n";
              //  //erase (tiles_str, )

              //  try_dup (ref, read, comstr, 
              //           f1, f2,
              //           g_hs, g_anchor, 
              //           dup_tiles_left_str, dup_tiles_rght_str,
              //           dup_str, dup_end, 
              //           thd_err_rate, 
              //           thd_tile_size,
              //           );
              //  reform_tiles(ref, read, comstr,
              //               dup_tiles_left_str, dup_tiles_left_end,
              //               dup_clips, dup_sp_tiles,
              //               dup_str, dup_end, g_map_rght,
              //               thd_cord_gap, thd_tile_size, thd_err_rate);
              //  reform_tiles(ref, read, comstr,
              //               dup_tiles_rght_str, dup_tiles_rght_end,
              //               dup_clips, dup_sp_tiles,
              //               dup_str, dup_end, g_map_left,
              //               thd_cord_gap, thd_tile_size, thd_err_rate); 
              //  g_print_tiles_(dup_tiles_rght_end, "dpts34");


              //  insert(tiles_str, getClipEnd(sp_tiles, i) + new_dups_count, dup_tiles_left_str);
              //  insert(tiles_end, getClipEnd(sp_tiles, i) + new_dups_count, dup_tiles_left_end);
              //  insert(tiles_str, getClipEnd(sp_tiles, i) + new_dups_count + length(dup_tiles_left_str), 
              //      dup_tiles_rght_str);
              //  insert(tiles_end, getClipEnd(sp_tiles, i) + new_dups_count + length(dup_tiles_left_end), 
              //      dup_tiles_rght_end);
              //  new_dups_count += length(dup_tiles_left_str) + length(dup_tiles_rght_str); 
                
            }
        }
    }
    else if (x1 > x2 && y1 < y2)
    {

    }
    else if (x1 < x2 && y1 > y2)
    {

    }
    else if (x1 > x2 && y1 > y2)
    {

    }
    //if (direction == 0)
    //{
        insertValue(tiles_str, 0, gap_str);
        insertValue(tiles_end, 0, shift_tile(gap_str, gap_parms.thd_tile_size, gap_parms.thd_tile_size));
        appendValue(tiles_str, shift_tile(gap_end, -gap_parms.thd_tile_size, -gap_parms.thd_tile_size));
        appendValue(tiles_end, gap_end);
        for (uint i = 1; i < length(tiles_str); i++)
        {
            int64_t dx = get_tile_x(tiles_str[i]) - get_tile_x(tiles_end[i - 1]);
            int64_t dy = get_tile_y(tiles_str[i]) - get_tile_y(tiles_end[i - 1]);
            if(get_tile_strand(tiles_str[i] ^ tiles_str[i - 1]))
            {

            }
            else
            {
                if (dx > 100 && dy > 100)
                {
                    String<uint64_t> tiles_str1;
                    String<uint64_t> tiles_end1;
                    String<uint64_t> sp_tiles_inv;
                    uint64_t t_gap_str = tiles_str[i - 1];
                    uint64_t t_gap_end = tiles_str[i];
                    /*
                    mapInterval(seqs[get_tile_id(gap_str)], read, comstr, tiles_str1, f1, f2,
                        t_gap_str, t_gap_end, LLMIN, LLMAX, t_direction, gap_parms);  
                    reform_tiles(seqs[get_tile_id(gap_str)], read, comstr, tiles_str1, tiles_end1, sp_tiles_inv, 
                        t_gap_str, t_gap_end, t_direction, gap_parms);
                     */   
                    mapGeneric (seqs, read, comstr, f1, f2, tiles_str1, tiles_end1, 
                        t_gap_str, t_gap_end, gap_parms);
                                //g_print_tiles_(tiles_end1, "mgs11");
                    if (!empty(tiles_str1))
                    {
                        erase(tiles_str1, 0); //inserted by reform_tiles
                        erase(tiles_end1, 0);
                        eraseBack(tiles_str1);
                        eraseBack(tiles_end1);
                        if (!empty(tiles_str1))
                        {
                            remove_tile_sgn(back(tiles_str1));
                            remove_tile_sgn(back(tiles_end1));
     //                       erase(tiles_str, i - 1);
     //                       erase(tiles_end, i - 1);
                            insert(tiles_str, i, tiles_str1);
                            insert(tiles_end, i, tiles_end1);
    //                        erase(tiles_str, i + length(tiles_str1) - 1);
    //                        erase(tiles_end, i + length(tiles_str1) - 1);
                        }
                        i += length(tiles_str1);
                    }
                    //dout << "mg2" << get_cord_y(gap_str) << get_cord_y(gap_end) << i << length(tiles_str) << dx << dy << gap_parms.read_len << direction << length(tiles_str1) << "\n";
                    //g_print_tiles_(tiles_str1, "gpt1");
                }
            }
        }
        for (int i = 1; i < (int)length(tiles_str) - 1; i++)
        {
            tiles_str[i - 1] = tiles_str[i];
            tiles_end[i - 1] = tiles_end[i];
        }
        resize(tiles_str, length(tiles_str) - 2);
        resize(tiles_end, length(tiles_end) - 2);
        //g_print_tiles_(tiles_end, "mgs3");
    //}
    return 0;
}
/**
 * Re-map gaps in cords.
 * Gaps at the front or end of the cords are also remapped.
 * !!NOTE::Each block in the cords are required to be sorted according to the y value of the first cord, such that isBlocksLinkable can work appropriately.
 */
int mapGaps(StringSet<String<Dna5> > & seqs, 
            String<Dna5> & read, 
            String<Dna5> & comstr,
            String<uint64_t> & cords_str, 
            String<uint64_t> & cords_end, 
            String<uint64_t> & clips, // string for clips cords
            String<UPair> & apx_gaps,
            StringSet<FeaturesDynamic > & f1,
            StringSet<FeaturesDynamic > & f2,
            GapParms & gap_parms)
{
    CmpInt64 g_cmpll;
    if (length(cords_str) <= 1)
    {
        return 0;
    }
    String <uint64_t> tiles_str;
    String <uint64_t> tiles_end;
    String<UPair> str_ends;
    String<UPair> str_ends_p;

    int64_t shift_x;
    int64_t shift_y;

    int thd_max_segs_num = 1000; //max segs num allowed in each gap, gaps > this will abort all tiles

    uint64_t thd_max_extend = 2000;  //Important::tune::whether process in gap or pmpfiner.h
    uint64_t thd_max_extend_x = thd_max_extend;
    uint64_t thd_max_extend_y = thd_max_extend; //>large gaps are supposed to be handled during the apx part.
    int64_t thd_max_gap = 3000;
    int64_t thd_dxy_min = 80;
    int64_t thd_extend_xy = 3; //extend at first or last cord of seqs
    int64_t thd_max_chain_distance = thd_max_extend;
    int64_t block_size = gap_parms.thd_tile_size; 
    int64_t thd_cord_size = gap_parms.thd_tile_size; 
    int64_t thd_cord_remap = 100;
    int64_t thd_cord_gap = gap_parms.thd_gap_len_min + block_size;

    unused(thd_cord_remap);
    unused(thd_max_chain_distance);

    clear(apx_gaps);
    gather_blocks_ (cords_str, str_ends, str_ends_p, 1, length(cords_str), length(read), thd_cord_gap, thd_cord_size, 0, &is_cord_block_end, &set_cord_end);
    gather_gaps_y_ (cords_str, str_ends, apx_gaps, length(read), thd_cord_gap);

    //NOTE cords_str[0] is the head, regular cords_str starts from 1
    for (unsigned i = 1; i < length(cords_str); i++)
    {
        unsigned sid = get_cord_id(cords_str[i]);
        int direction;
        gap_parms.read_len = length(read);
        gap_parms.ref_len = length(seqs[sid]);

        if (_DefaultCord.isBlockEnd(cords_str[i - 1]))  //clip first cord
        {
            shift_x = std::min (int64_t(length(seqs[sid]) - 1 - get_cord_x(cords_str[i])),  block_size);
            shift_y = std::min (int64_t(length(read) - 1 - get_cord_y(cords_str[i])), block_size);
            uint64_t gap_end = shift_cord(cords_str[i], shift_x, shift_y);
            uint64_t gap_str;
            if ((int64_t)get_cord_y(gap_end) > thd_cord_gap)
            {            
                shift_x = std::min(thd_max_extend, get_cord_x(gap_end));
                shift_y = std::min(thd_max_extend, get_cord_y(gap_end));
                shift_x = std::min(shift_x, shift_y * thd_extend_xy);
                uint64_t infi_cord = shift_cord(gap_end, -shift_x, -shift_y); 
                direction = g_map_left;
                gap_str = infi_cord;
                _DefaultHit.unsetBlockEnd(gap_str);
                _DefaultHit.unsetBlockEnd(gap_end);

                int max_gap_overlap_y = _getMaxGapsyOverlap(apx_gaps, gap_str, gap_end);
                if (max_gap_overlap_y > thd_cord_gap)
                {
                    mapGap_ (seqs, read, comstr, gap_str, gap_end, f1, f2, tiles_str, tiles_end, clips, direction, thd_dxy_min, gap_parms);
                    //g_print_tiles_(tiles_end, "mg1");
                    insert_tiles2Cords_(cords_str, cords_end, i, tiles_str, tiles_end, direction, thd_cord_size, thd_max_segs_num);
                    //<<debug
                    //g_print_tiles_(tiles_end, "mg2");
                    //>>debug
                }
                //_DefaultHit.setBlockEnd(cords_str[i - 1 + length(tiles_str)]);
                //_DefaultHit.setBlockEnd(cords_end[i - 1 + length(tiles_str)]);
                /*
                else
                {
                    uint64_t tile_str = cords_str[i];
                    uint64_t tile_end = shift_cord(cords_str[i], block_size, block_size);

                    reform_tile_ (seqs[get_cord_id(cords_str[i])], read, comstr,
                                  tile_str, tile_end, g_sv_l, thd_tile_size);
                    _updateCordsStrEndValue(cords_str, cords_end, i, tile_str, tile_end, block_size);
                    
                    int64_t anchor_base = g_hs_Cord2StrAnchor(gap_end);
                    int64_t anchor_lower = anchor_base - 100;
                    int64_t anchor_upper = anchor_base + 100;
                    print_cord(gap_str, "sstr1");
                    print_cord(gap_end, "sstr2");
                    mapExtend_ (seqs, read, comstr,  
                                f1, f2, tiles_str, tiles_end, clips, 
                                gap_str, gap_end, direction,
                                thd_cord_gap, thd_tile_size, thd_cord_remap, thd_err_rate, thd_dxy_min, gap_parms);
                    insert_tiles2Cords_(cords_str, cords_end, i, tiles_str, tiles_end, direction, thd_cord_size, thd_max_segs_num);
                }
                */
            }
        }
        else if (!isCordsConsecutive_(cords_str[i - 1], cords_str[i], thd_cord_gap))
        {
            //continue;
            g_cmpll.min(shift_x, block_size) << int64_t(length(seqs[sid]) - 1 - get_cord_x(cords_str[i]));
            g_cmpll.min(shift_y, block_size) << int64_t(length(read) - 1 - get_cord_y(cords_str[i]));
            uint64_t gap_end = shift_cord(cords_str[i], shift_x, shift_y);
            uint64_t gap_str = cords_str[i - 1]; 
            int64_t dx_tmp = get_cord_x(gap_end) - get_cord_x(gap_str);
            if (std::abs(int64_t(dx_tmp)) < thd_max_gap)
            {
            //dout << "g1" << get_cord_y(gap_str) << get_cord_y(gap_end) << get_cord_x(gap_str) << get_cord_x(gap_end) << dx_tmp << i - 1 << i << "\n";
                direction = g_map_closed;
                _DefaultHit.unsetBlockEnd(gap_str);
                _DefaultHit.unsetBlockEnd(gap_end);
                
                mapGap_(seqs, read, comstr, gap_str, gap_end, f1, f2, tiles_str, tiles_end,
                        clips, direction, thd_dxy_min, gap_parms);
                insert_tiles2Cords_(cords_str, cords_end, i, tiles_str, tiles_end, direction, thd_cord_size, thd_max_segs_num);
            }
        }
        if (_DefaultHit.isBlockEnd(cords_str[i]))  ///right clip end cord
        {
            //continue;
            uint64_t gap_str = cords_str[i];
            if (int64_t(length(read) - 1 - get_cord_y(gap_str)) > thd_cord_gap)
            {
                shift_x = std::min(thd_max_extend_x, length(seqs[sid]) - get_cord_x(gap_str) - 1);
                shift_y = std::min(thd_max_extend_y, length(read) - get_cord_y(gap_str) - 1);
                shift_x = std::min(shift_x, shift_y * thd_extend_xy);
                uint64_t infi_cord = _DefaultCord.shift(gap_str, shift_x, shift_y);
                uint64_t gap_end = infi_cord;
                direction = g_map_rght;
                _DefaultHit.unsetBlockEnd(gap_str);
                _DefaultHit.unsetBlockEnd(gap_end);
                int max_gap_overlap_y = _getMaxGapsyOverlap(apx_gaps, gap_str, gap_end);
                if (max_gap_overlap_y > thd_cord_gap)
                {
                    mapGap_ (seqs, read, comstr, gap_str, gap_end, f1, f2, tiles_str, tiles_end, clips, direction, thd_dxy_min, gap_parms);
                    insert_tiles2Cords_(cords_str, cords_end, i, tiles_str, tiles_end, direction, thd_cord_size, thd_max_segs_num);
                    //dout << "mgr" << is_cord_block_end(cords_str[i + length()]) << "\n";
                }
                /*else
                {
                    uint64_t tile_str = cords_str[i];
                    uint64_t tile_end = shift_cord(cords_str[i], block_size, block_size);
                    reform_tile_ (seqs[get_cord_id(cords_str[i])], read, comstr, tile_str, tile_end, g_sv_r,thd_tile_size);
                    _updateCordsStrEndValue(cords_str, cords_end, i, tile_str, tile_end, block_size);

                    int64_t anchor_base = g_hs_Cord2StrAnchor(gap_str);
                    int64_t anchor_lower = anchor_base - 100;
                    int64_t anchor_upper = anchor_base + 100; 
                    print_cord(gap_str, "eedn1");
                    print_cord(gap_end, "eedn2");
                    print_cord(cords_str[i], "eedn3");
                    mapExtend_ (seqs, read, comstr,  
                                f1, f2, tiles_str, tiles_end, clips, 
                                gap_str, gap_end, direction,
                                thd_cord_gap, thd_tile_size, thd_cord_remap, thd_err_rate, thd_dxy_min, gap_parms);
                    insert_tiles2Cords_(cords_str, cords_end, i, tiles_str, tiles_end, direction, thd_cord_size, thd_max_segs_num);
                }
                */
            }
        }
    }
    //mapGapGlobal();
    return 0;
}

/*=====  End of Index free Map and clip  ======*/



