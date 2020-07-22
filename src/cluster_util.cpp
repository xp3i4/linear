#include <iostream>
#include "base.h"
#include "cords.h"
#include "cluster_util.h"

using namespace seqan;
using std::cout;
using std::endl;

/*__________________________________________________
  ---------- @s::Generic funcs  ----------*/
ChainScoreParms::ChainScoreParms()
{
    mean_d = 1000; 
    var_d = 1000;
    chn_block_strand = 0;
};

int ChainsRecord::isLeaf(){return f_leaf;};

//Chainning Score metric wrapper: including a score function with corresponding parms.
ChainScoreMetric::ChainScoreMetric(){};
ChainScoreMetric::ChainScoreMetric(int min_chain_len, int abort_socre, 
        int(*scoreFunc)(uint64_t const &, uint64_t const &, ChainScoreParms &)) 
        : thd_min_chain_len(min_chain_len), thd_abort_score(abort_socre), getScore(scoreFunc), getScore2(NULL)
        {};

ChainScoreMetric::ChainScoreMetric(int min_chain_len, int abort_score, int (*scoreFunc)(uint64_t const &, 
    uint64_t const &, uint64_t const &, uint64_t const &, uint64_t const & read_len, ChainScoreParms &)) :
    thd_min_chain_len(min_chain_len),
    thd_abort_score(abort_score), getScore(NULL), getScore2(scoreFunc)
    {};    
int ChainScoreMetric::getMinChainLen()
{
    return thd_min_chain_len;
}
int ChainScoreMetric::getAbortScore()
{
    return thd_abort_score;
}

int const chain_end_score = 0;
int const chain_end = -1;

/*
 * Warn:: anchors are required to be sorted by x in descending order.
   Don't manipulate the anchors by direct comparation, insertion in this function !
   The @anchors is an abstract class (not defined), which can be cords, gaps, or hit.
   It can only be handled by the corresponding @scoreFunc, which defines the @anchors.
 * chain anchors within [@it_str, @it_end)
 */ 
int getBestChains(String<uint64_t>     & anchors, //todo:: anchor1 anchor2 of different strand not finished 
                  String<ChainsRecord> & chains,
                  uint it_str, uint it_end, uint64_t thd_chain_depth, uint64_t thd_chain_dx_depth,
                  ChainScoreMetric & chn_metric,
                  uint64_t (*get_anchor_x)(uint64_t))
{
    if (empty(anchors))
    {
        return 0;
    }
    int new_score = 0;
    int new_max_score = 0;
    int max_j = 0;
    int const chain_end_score = 0;
    int const chain_end = -1;

    chains[0].score = chain_end_score;  
    chains[0].len = 1;
    chains[0].p2anchor = chain_end;
    //for (int i = 0; i < length(anchors); i++) 
    for (int i = it_str; i < it_end; i++) 
    {
        int j_str = std::max (0, i -  (int)thd_chain_depth);
        max_j = i;
        new_max_score = -1;
        //for (int j = j_str; j < i; j++)
        for (int j = i - 1; j>=0 && (j >=j_str || get_anchor_x(anchors[j]) - get_anchor_x(anchors[i]) < thd_chain_dx_depth); j--)
        {
            new_score = chn_metric.getScore(anchors[j], anchors[i], chn_metric.chn_score_parms);
            if (new_score > 0 && new_score + chains[j].score >= new_max_score)
            {
                max_j = j;
                new_max_score = new_score + chains[j].score;
            }
        }
        if (new_max_score > 0)
        {
            chains[i].p2anchor = max_j;
            chains[i].score = new_max_score ;
            chains[i].len = chains[max_j].len + 1;
            chains[i].score2 = new_max_score;
            chains[i].root_ptr = chains[max_j].root_ptr;
            chains[i].f_leaf = 1;
            chains[max_j].f_leaf = 0;
        }
        else
        {
            chains[i].p2anchor = chain_end;
            chains[i].score = chain_end_score;
            chains[i].len = 1;
            chains[i].score2 = chain_end_score;
            chains[i].root_ptr = i;
            chains[i].f_leaf = 1;
        }
    }
    return 0;
}

/* Back trace function to retrive @chains of @elements from the @chains_record
 * For any element in the @elements, it will be chained at most once in the most likely(highest score) chain
 */
template <class ChainElementType>
int traceBackChains0(String<ChainElementType> & elements,  StringSet<String<ChainElementType> > & chains, String<ChainsRecord> & chain_records, String<int> & chains_score, int _chain_min_len, int _chain_abort_score, int bestn)
{
    String<ChainElementType> chain;
    String<int> chain_score;
    int delete_score = -1000;
    int search_times = 8;
    int num_chains = 0;
    for (int i = 0; i < search_times && num_chains++ <= bestn; i++) 
    {
        bool f_done = true;
        int max_2nd_score = -1;
        int max_score = -1;
        int max_str = chain_end;
        int max_len = 0;
        for (int j = 0; j < length(chain_records); j++)
        {
            if (chain_records[j].score > max_score)
            {
                max_2nd_score = max_score;
                max_str   = j;
                max_score= chain_records[j].score;
                max_len   = chain_records[j].len;
                f_done = false;
            }
        }
        if (f_done || max_score == 0)
        {
            break;
        }
        if (max_len > _chain_min_len && max_score / (max_len - 1) > _chain_abort_score) //max_len is the number of anchors, ..-1 is the number of connection(interval) between anchors
        {
            for (int j = max_str; j != chain_end; j = chain_records[j].p2anchor)
            {
                if (chain_records[j].score != delete_score)
                {
                    appendValue (chain, elements[j]);
                    appendValue (chain_score, chain_records[j].score2);
                    chain_records[j].score = delete_score; 
                }
                else
                {
                    int infix_chain_score = chain_records[j].score2; //infix of the chain has been seleted out, so update the score of the suffix chain.

                    if (max_score - infix_chain_score < max_2nd_score)
                    {
                        for (int k = max_str; k != j; k = chain_records[k].p2anchor)  
                        {
                            chain_records[k].score = chain_records[k].score2 - infix_chain_score;
                        } 
                        clear (chain); 
                        clear (chain_score);
                    }

                    break;
                }
            }
            if (!empty(chain))
            {
                appendValue(chains, chain);
                append(chains_score, chain_score);
                clear(chain);
                clear(chain_score);
            }
        } 
        if (max_str != chain_end)
        {
            chain_records[max_str].score = delete_score; 
        }
    }  
    return 0;
}

/* Back trace function to retrive @chains of @elements from the @chains_record
 * For any element in the @elements, it will be chained at most once in the most likely(highest score) chain
 */
template <class ChainElementType>
int traceBackChains(String<ChainElementType> & elements,  StringSet<String<ChainElementType> > & chains, String<ChainsRecord> & chain_records, String<int> & chains_score, int _chain_min_len, int _chain_abort_score, int bestn)
{
    String<ChainElementType> chain;
    String<int> chain_score;
    int delete_score = -1000;
    int search_times = 8;
    int num_chains = 0;
    String<int> new_leaves;
    StringSet<String<int> > leaves;
    resize (new_leaves, 4);

    for (int j = 0; j < length(chain_records); j++)
    {
        if (chain_records[j].isLeaf()) //create leaves list for each tree
        {
            int f_new = 1;
            for (int k = 0; k < length(leaves); k++)
            {
                if (leaves[k][0] == chain_records[j].root_ptr)
                {
                    appendValue(leaves[k], j);
                    if (chain_records[j].score > leaves[k][1])
                    {
                        leaves[k][1] = chain_records[j].score;
                        leaves[k][2] = chain_records[j].len;
                        leaves[k][3] = j;
                    }
                    f_new = 0;
                }
            }
            if (f_new)
            {
                new_leaves[0] = chain_records[j].root_ptr;
                new_leaves[1] = chain_records[j].score;
                new_leaves[2] = chain_records[j].len;
                new_leaves[3] = j; //leaf of max score
                appendValue(leaves, new_leaves);
            }
        }
    }
    
    String<std::pair<int,int> > tree_score_ranks;
    resize(tree_score_ranks, length (leaves));
    for (int i = 0; i < length(leaves); i++)
    {
        tree_score_ranks[i] = std::pair<int, int> (i, leaves[i][1]);
    }
    std::sort (begin(tree_score_ranks), end(tree_score_ranks), 
        [](std::pair<int, int> & a, std::pair<int, int> & b){return a.second > b.second;});
    for (int i = 0; i < std::min(bestn, int(length(tree_score_ranks))); i++) 
    {
        int max_score = leaves[tree_score_ranks[i].first][1];
        int max_len = leaves[tree_score_ranks[i].first][2];
        int max_str = leaves[tree_score_ranks[i].first][3];
        int mean_score = max_len > 1 ? max_score / (max_len - 1) : _chain_abort_score + 1;
        if (max_len > _chain_min_len && mean_score > _chain_abort_score) //max_len is the number of anchors, ..-1 is the number of connection(interval) between anchors
        {
            for (int j = max_str; j != chain_end; j = chain_records[j].p2anchor)
            {
                appendValue (chain, elements[j]);
                appendValue (chain_score, chain_records[j].score2);
            }
            if (!empty(chain))
            {
                appendValue(chains, chain);
                append(chains_score, chain_score);
                clear(chain);
                clear(chain_score);
            }
        } 
    }  
    return 0;
}

int getApxChainScore(uint64_t const & anchor1, uint64_t const & anchor2, ChainScoreParms & chn_sc_parms)
{
    int64_t dy = get_cord_y(anchor1) - get_cord_y(anchor2);
    if (dy < 10)
    {
        //dy < 0 : y should in descending order
        //0 <= dy < 10 : too close anchors are excluded;
        return -10000;
    }
    int64_t thd_min_dy = 50;
    //int64_t dx = dy + int64_t(_DefaultHit.getAnchor(anchor2) - _DefaultHit.getAnchor(anchor1));
    int64_t dx = getAnchorX(anchor1) -  getAnchorX(anchor2);
    int64_t da = std::abs(dx - dy);
    int64_t derr =  (100 * da) / std::max({std::abs(dy), std::abs(dx), thd_min_dy}); // 1/100 = 0.01
    
    //d_err
    /*
    if (derr < 10)         {derr = 0;}
    else if (derr < 25)    {derr = 10 + 2 * derr ;}
    else if (derr < 100)   {d_err = derr * derr / 10 + 40;}
    else                    {derr = 10000;}
    */
    //d_err
    int score_derr;
    if (derr < 5)
    {
        score_derr = 4 * derr;
    }
    else if (derr < 10)
    {
        score_derr = 6 * derr - 10;
    }
    else if (derr < 100) 
    {
        score_derr =  derr * derr - 5 * derr;
    }
    else
    {
        return -1000;
    }
    //d_y
    int score_dy;
    dy /= 15;
    if (dy < 150)           {score_dy = dy / 5;}
    else if (dy < 100)      {score_dy = dy - 30;}
    else if (dy < 10000)    {score_dy = dy * dy / 200 + 20;}
    else                    {score_dy = 10000;}
    if (da < 10)
    {
        return 100 - score_dy;
    }
    else
    {
        return 100 - score_dy - score_derr ;    
    }
}

int chainAnchorsBase(String<uint64_t> & anchors, StringSet<String<uint64_t> > & anchors_chains, 
    String<int> & anchors_chains_score, uint it_str, uint it_end,  uint thd_chain_depth, uint64_t thd_chain_dx_depth, 
    int thd_best_n, ChainScoreMetric & chn_metric, uint64_t (*get_anchor_x) (uint64_t))
{
    if (length(anchors) < 2){
        return 0;
    }
    String<ChainsRecord> chain_records;
    resize (chain_records, length(anchors));
    getBestChains (anchors, chain_records, it_str, it_end, thd_chain_depth, thd_chain_dx_depth, chn_metric, get_anchor_x);
    traceBackChains(anchors, anchors_chains, chain_records, anchors_chains_score, chn_metric.getMinChainLen(), chn_metric.getAbortScore(), thd_best_n);
    return 0;
}

int chainAnchorsHits(String<uint64_t> & anchors, String<uint64_t> & hits, String<int> & hits_chains_score)
{
    int thd_min_chain_len = 1;
    int thd_drop_score = 45; //<<TODO, change the score!
    uint thd_chain_depth = 20;
    ChainScoreMetric chn_score(thd_min_chain_len, thd_drop_score, &getApxChainScore);
    StringSet<String<uint64_t> > anchors_chains;
    std::sort(begin(anchors), end(anchors), 
        [](uint64_t & a, uint64_t & b){return getAnchorX(a) > getAnchorX(b);});
    int thd_best_n = 5;
    chainAnchorsBase(anchors, anchors_chains, hits_chains_score, 0, length(anchors), thd_chain_depth, 0, thd_best_n, chn_score, &getAnchorX);
    //additoinal filter and convert to hits
    for (int i = 0; i < length(anchors_chains); i++)
    {
        for (int j = 0; j <length(anchors_chains[i]); j++)
        {
            appendValue(hits, _DefaultCord.hit2Cord_dstr(anchors_chains[i][j]));
        }
        _DefaultHit.setBlockEnd(back(hits));
    } 
    return 0;
}

/*
 * This is the copy version of getBestChains() for ChainBlock using the same algorithm.
 * Note::So synchronsize with the getBesctChains().
 * For efficiency when profiling, just duplicate the func instead of using template or virtual function
 */
int getBestChains2(String<uint64_t> & hits,
                   String<UPair> & str_ends_p,
                   String<int>   & str_ends_p_score,
                   String<ChainsRecord> & chain_records,
                   uint64_t read_len,
                   ChainScoreMetric & chn_metric)
{
    int thd_chain_depth = 20;
    int new_score = 0;
    int new_max_score = 0;
    int max_j = 0;
    int const chain_end_score = 0;
    int const chain_end = -1;

    chain_records[0].score = str_ends_p_score[0];  
    chain_records[0].len = str_ends_p[0].second - str_ends_p[0].first;
    chain_records[0].p2anchor = chain_end;
    for (int i = 0; i < length(str_ends_p); i++) 
    {
        int j_str = std::max (0, i - thd_chain_depth);
        max_j = i;
        new_max_score = -1;
        for (int j = j_str; j < i; j++)
        {
            new_score = chn_metric.getScore2
                            (hits[str_ends_p[j].first], hits[str_ends_p[j].second - 1],
                             hits[str_ends_p[i].first], hits[str_ends_p[i].second - 1],
                                   read_len, chn_metric.chn_score_parms);

            if (new_score > 0 && new_score + chain_records[j].score + str_ends_p_score[i] >= new_max_score)
            {
                max_j = j;
                new_max_score = new_score + chain_records[j].score + str_ends_p_score[i];
            }
        }
        if (new_max_score > 0)
        {
            chain_records[i].p2anchor = max_j;
            chain_records[i].score = new_max_score ;
            chain_records[i].len = str_ends_p[i].second - str_ends_p[i].first + chain_records[max_j].len;
            chain_records[i].score2 = chain_records[i].score;
            chain_records[i].root_ptr = chain_records[max_j].root_ptr;
            chain_records[i].f_leaf = 1;
            chain_records[max_j].f_leaf = 0;
        }
        else
        {
            chain_records[i].p2anchor = chain_end;
            chain_records[i].score = str_ends_p_score[i];
            chain_records[i].len = str_ends_p[i].second - str_ends_p[i].first;
            chain_records[i].score2 = chain_records[i].score;
            chain_records[i].root_ptr = i;
            chain_records[i].f_leaf = 1;
        }
    }
    return 0;
}

/*
 * @scoreFunc is the score function to create chains of blocks of records of @records
 * @records is supposed to be in the format of cords
 * @chains[0] is the best chain
 */
int chainBlocksBase(StringSet<String<UPair> > & chains, String<uint64_t> & records, String<UPair> & str_ends_p, String<int> & str_ends_p_score, uint64_t read_len, ChainScoreMetric & chn_metric, int thd_best_n, int f_sort = 1)
{
    if (length(str_ends_p) < 2) {
        return 0;
    }
    String<ChainsRecord> chain_records;
    String<int> chains_score;
    String <unsigned> ptr;
    String<UPair> str_ends_p_tmp;
    String<int> str_ends_p_score_tmp;
    //sort str_ends_p and str_ends_p_score in denscending of corresponding x of records[str_end_p]; ptr is tmp pointer array
    for (int i = 0; i < length(str_ends_p); i++)
    {
        appendValue(ptr, i);
    }
    if (f_sort)
    {
        std::sort (begin(ptr), end(ptr), [& records, & str_ends_p](unsigned & a, unsigned & b){
            return _DefaultCord.getCordX(records[str_ends_p[a].first]) > _DefaultCord.getCordX(records[str_ends_p[b].first]);
        });
    }
    resize(str_ends_p_tmp, length(str_ends_p));
    resize(str_ends_p_score_tmp, length(str_ends_p_score));
    for (int i = 0; i < length(str_ends_p); i++)
    {
        str_ends_p_tmp[i] = str_ends_p[ptr[i]];
        str_ends_p_score_tmp[i] = str_ends_p_score[ptr[i]];
    }

    resize (chain_records, length(str_ends_p_tmp));
    getBestChains2(records, str_ends_p_tmp, str_ends_p_score_tmp, chain_records, read_len, chn_metric);
    traceBackChains(str_ends_p_tmp, chains, chain_records, chains_score, chn_metric.getMinChainLen(), chn_metric.getAbortScore(), thd_best_n);
    return 0;
}

/*
 * Warn:: x
 * score of chain block [@cord11, @cord12) and block [@cord21, @cord22)
 * @cord21 chain to @cord11
 * x22 < x11 are required
 * @cord*1 and @cord*2 are required to have the same strand
 */
int getApxChainScore2(uint64_t const & cord11, uint64_t const & cord12, uint64_t const & cord21, uint64_t const & cord22, uint64_t const & read_len, ChainScoreParms & chn_sc_parms)
{
    int64_t thd_max_d = 20000;
    int64_t thd_indel_trigger = 100;
    int64_t thd_indel_op = 30; //indel open penalty
    int64_t dy = get_cord_y(cord11) - get_cord_y(cord22);
    int64_t dx = get_cord_x(cord11) - get_cord_x(cord22);
    if (dx < 0 || dy < 0 || get_cord_strand(cord11 ^ cord22) || dx > thd_max_d || dy > thd_max_d)
    {
        return INT_MIN;
    }
    int64_t thd_min_dy = 100;
    int64_t da = std::abs(int64_t(dx - dy));
    int64_t derr =  (100 * da) / std::max({std::abs(dy), thd_min_dy, std::abs(dx)}); // 1/100 = 0.01
    int score_derr;
    if (da > thd_indel_trigger || derr > 50)
    {
        if (dx < dy) //ins
        {
            return 100 - thd_indel_op - dy / 1000 - dx / 100;
        }
        else //del
        {
            return 100 - thd_indel_op - dy /100 - dx / 1000;
        }

    }
    else
    {
        return 100 - dy / 95;
    }
    /*
    if      (d_err < 10)    {d_err = 0;}                 //10%
    else if (d_err < 25)    {d_err = 10 + d_err;}        //25%
    else if (d_err < 100)   {d_err = 35 + d_err / 2;}
    else                    {d_err = 10000;}
    */
    //dy /= 150;    

    (void)read_len;

    //return 100 - dy - score_derr ;    
}

int _filterBlocksHits(StringSet<String<UPair> > & chains, String<uint64_t> & hits, uint64_t read_len)
{
    if (empty(chains))
    {
        return 0;
    }
    //step 2 filter major chain and remove poorly chained hits blocks
    //chains[0] is the major chain
    String<UPair> best_chain;
    String<uint64_t> hits_tmp;
    uint64_t len_current = 0;
    resize(best_chain, length(chains[0]));
    for (uint i = 0; i < length(chains[0]); i++)
    {
        for (uint j = chains[0][i].first; j < chains[0][i].second; j++)
        {
            appendValue(hits_tmp, hits[j]);
            _DefaultHit.unsetBlockEnd(back(hits_tmp));
        }
        len_current += chains[0][i].second - chains[0][i].first;
        best_chain[i] = chains[0][i];
    }
    _DefaultHit.setBlockEnd(back(hits_tmp));
    //process the chains left
    float thd_major_bound = 0.8 * len_current; // len > this * first major len is regarded as optional major chain
    uint thd_major_limit = 2;
    uint major_n = 1;
    uint64_t thd_x_max_delta = read_len * 2; //max distance allowed any x to the x of the major chain
    bool f_append = false;
    for (uint i = 1; i < length(chains); i++)
    {
        len_current = 0;
        f_append = false;
        for (uint j = 0; j < length(chains[i]); j++)
        {
            len_current += chains[i][j].second - chains[i][j].first; 
        }
        if (major_n < thd_major_limit && len_current > thd_major_bound) // the 2nd,3th.. optional major chain
        {
            f_append = true;
            ++major_n;
        }
        else //check if can append(co-exists) to the 1st major chain (inversions etc..) 
        {
            f_append = true;
            for (int j = 0; j < length(chains[i]) && f_append; j++)
            {
                for (int k = 0; k < length(best_chain) && f_append; k++) 
                {
                    uint64_t str_major = hits[best_chain[k].first];
                    uint64_t end_major = hits[best_chain[k].second - 1];
                    uint64_t str_current = hits[chains[i][j].first];
                    uint64_t end_current = hits[chains[i][j].second - 1];
                    int64_t dx_lower = int64_t(get_cord_x(str_major) - get_cord_x(str_current));
                    int64_t dx_upper = int64_t(get_cord_x(end_current) - get_cord_x(end_major));
                    f_append = dx_lower <= thd_x_max_delta && dx_upper < thd_x_max_delta  &&
                               !_isCordyOverLap(str_major, end_major, str_current, end_current, read_len);
                }
            }
            for (int j = 0; j < length(chains[i]) && f_append; j++) //if append then extend the best_chain
            {
                append(best_chain, chains[i]);
            }        
        }
        if (f_append)
        {
            for (int j = 0; j < length(chains[i]); j++)
            {
                for (int k = chains[i][j].first; k < chains[i][j].second; k++)
                {
                    appendValue(hits_tmp, hits[k]);
                    _DefaultHit.unsetBlockEnd(back(hits_tmp));
                }
            }
            _DefaultHit.setBlockEnd(back(hits_tmp));
        }
            _DefaultHit.setBlockEnd(back(hits_tmp));
    }
    hits = hits_tmp;  
    return 0;
}

int chainBlocksHits(String<uint64_t> & hits, String<UPair> & str_ends_p, String<int> & str_ends_p_score, uint64_t read_len)
{
    int thd_min_chain_len = 1;
    int thd_drop_score = 0;
    int thd_best_n = 3;
    StringSet<String<UPair> > hits_chains;
    ChainScoreMetric chn_score(thd_min_chain_len, thd_drop_score, &getApxChainScore2); //init the score as the length of the blocks
    chainBlocksBase(hits_chains, hits, str_ends_p, str_ends_p_score, read_len, chn_score, thd_best_n);
    _filterBlocksHits(hits_chains, hits, read_len);
    return 0;
}

/*
 * Comment of getForwardChainDxDy and getApxChainScore3
 * Loosely chaining of cords blocks of different varints
 * cord blocks chaining  [+15, +16), [-80, -81), [-82, -83),[+21, +22), wher + - as forward and reversed 
 *          o1     o2     o3
    + 15 | + 15 | + 21 s| + 21  y11 
      16 |   16 |   22  |   22  y12
    -----|------|-------|------
    - 80 | + 20 | + 20  | - 80  y21
      81 |   19 |   19 s|   81  y22
    -----|------|-------|------ 
    - 82 | + 18 | + 18  | - 82
      83 |   17 |   17 s|   83
    -----|------|-------|------
    + 21 | + 21 | + 15 s| + 15
      22 |   22 |   16  |   16
   o1: y to the forward strand
   o2:sort each block by start_y of on forwardstrand in descending order
   o3:o1*o3=1, y to original strand
   element marked 's' is the key to sort
   result in column o3 is passed to this function

 *[Warn::red] the fucntion requires cords to be pre-sorted by y rather than x, that is different from the other chaining functions.
 * Last stage of chaining on level cords already extended. 
 * score of chain cords block [@cord11, @cord12) and block [@cord21, @cord22)
 * @cord21 chain to @cord11
 * @cord*1 and @cord*2 is allowed to have the different strands in case of inversions
 * @dx is allowed to be < 0 in case of duplications 
 _____________________
 * STRATEGY::
   1.regular cords has the highest priority to be chainned
   2.variants signals including ins,dup,inv,del has the same priority to be chainned.
 * Drawbacks: in some cases dup can't be chained, e.g, 
   ---can't since 9 < 10
   block1, (10, 10) - (15, 15)
   block2, (9,  16) - (13, 20)
   ---can 
   block1, (10, 10) - (15, 15)
   block2, (11, 18) - (13, 20) 
*/
int getChainBlockDxDy(uint64_t const & cord11, uint64_t const & cord12, uint64_t const & cord21, uint64_t const & cord22, uint64_t const & read_len, int strand, int64_t & dx, int64_t & dy)
{
    uint64_t y1, y2; 
    if (get_cord_strand(cord11) != strand)
    {
        if (get_cord_strand(cord22) != strand)
        {
            dy = get_cord_y(cord21) - get_cord_y(cord12);
            dx = get_cord_x(cord21) - get_cord_x(cord12);
        }
        else
        {
            dy = read_len - get_cord_y(cord12) - 1 - get_cord_y(cord22);
            dx = get_cord_x(cord11) - get_cord_x(cord22);
        }
    }
    else
    {
        if (get_cord_strand(cord22) != strand)
        {
            dy = get_cord_y(cord11) - read_len + 1 + get_cord_y(cord21);
            dx = get_cord_x(cord11) - get_cord_x(cord22);
        }
        else
        {
            dy = get_cord_y(cord11) - get_cord_y(cord22);
            dx = get_cord_x(cord11) - get_cord_x(cord22);
        }
    }

    return get_cord_strand(cord11 ^ cord22);
}

//Warn red>sychronize getGapChainScore2 of same logic  when modifiy this function  
int getApxChainScore3(uint64_t const & cord11, uint64_t const & cord12, uint64_t const & cord21, uint64_t const & cord22, uint64_t const & read_len, ChainScoreParms & chn_sc_parms)
{
    int64_t thd_min_dy = -80;
    int64_t thd_min_dx = -80;
    int64_t dx, dy, da, d_err; 
    //int f_type = getForwardChainDxDy(cord11, cord12, cord21, cord22, read_len, dx, dy);
    int f_type = getChainBlockDxDy(cord11, cord12, cord21, cord22, read_len, chn_sc_parms.chn_block_strand, dx, dy);
    
    int64_t thd_max_dy = 3000; 
    int64_t thd_max_dx = 15000;
    int64_t thd_dup_trigger = -50;
    int64_t dx_ = std::abs(dx);
    int64_t dy_ = std::abs(dy);
    da = dx - dy;
    int score = 0;
    if (dy < thd_min_dy)// || (f_type == 0 && dy > thd_max_dy) || dx_ > thd_max_dx)
    {
        score = INT_MIN;
    }
    else
    {
        int64_t score_dy = dy_ > 2000 ? std::min(dy_ / 25 - 50, int64_t(70)): dy_ / 40;  
        int64_t score_dx = dx_ > 2000 ? std::min(dx_ / 25 - 50, int64_t(70)): dx_ / 40;  
        if (f_type == 1) //inv
        {
            if (dx > thd_min_dx)
            {
                score = 75 - score_dy; 
            }
        }
        else if (da < -std::max(dx_ / 4, int64_t(50))) 
        {
            if (dx > thd_dup_trigger) //ins
            {
                score = 80 - score_dx; // large dy is allowed
            }
            else //dup
            {
                //todo limit dx < read_len
                score = 80 - score_dy; //  dy of dup is suppoesd to be close enough
            }
        }
        else if (da > std::max(dy / 4, int64_t(50))) //del
        {
            score = 80 - score_dy;
        }
        else //normal 
        {
            score = 100 - score_dy;
        }
    }
    dout << "cs3" << dx << dy << chn_sc_parms.chn_block_strand << get_cord_y(cord11) << get_cord_y(cord22) << get_cord_x(cord11) << get_cord_x(cord22) << score << "\n";
    return score;
}

int _filterBlocksCords(StringSet<String<UPair> > & chains, String<uint64_t> & hits, uint64_t read_len, uint64_t thd_major_limit, 
        void (*unsetEndFunc)(uint64_t &), void(*setEndFunc)(uint64_t &), int f_header)
{
    if (empty(chains))
    {
        return 0;
    }
    //step 2 filter major chain and remove poorly chained hits blocks
    //chains[0] is the major chain
    String<UPair> best_chain;
    String<uint64_t> hits_tmp;
    uint64_t len_current = 0;
    resize(best_chain, length(chains[0]));
    if (f_header) //cords[0] is header, while tiles[0] is not
    {
        appendValue (hits_tmp, hits[0]);
    }   
    for (uint i = 0; i < length(chains[0]); i++)
    {
        for (uint j = chains[0][i].first; j < chains[0][i].second; j++)
        {
            appendValue(hits_tmp, hits[j]);
            unsetEndFunc(back(hits_tmp));
            //_DefaultHit.unsetBlockEnd(back(hits_tmp));
        }
        len_current += chains[0][i].second - chains[0][i].first;
        best_chain[i] = chains[0][i];
    }
    setEndFunc(back(hits_tmp));
    //_DefaultHit.setBlockEnd(back(hits_tmp));
    //process the chains left
    float thd_major_bound = 0.8 * len_current; // len > this * first major len is regarded as optional major chain
    uint major_n = 1;
    uint64_t thd_x_max_delta = read_len * 2; //max distance allowed any x to the x of the major chain
    bool f_append = false;
    for (uint i = 1; i < length(chains) && major_n < thd_major_limit; i++)
    {
        len_current = 0;
        f_append = false;
        for (uint j = 0; j < length(chains[i]); j++)
        {
            len_current += chains[i][j].second - chains[i][j].first; 
        }
        if (len_current > thd_major_bound) // the 2nd,3th.. optional major chain
        {
            f_append = true;
            ++major_n;
        }
        if (f_append)
        {
            for (int j = 0; j < length(chains[i]); j++)
            {
                for (int k = chains[i][j].first; k < chains[i][j].second; k++)
                {
                    appendValue(hits_tmp, hits[k]); 
                    unsetEndFunc(back(hits_tmp));
                }
            }
            setEndFunc(back(hits_tmp));
        }
    }
    hits = hits_tmp;  
    return 0;
}
/*
 * Chain blocks on one strand specified by @chn_score.chn_score_parms.chn_block_strand
 * y of blocks of different strands are converted to y' of the same strand and chained
 */
int chainBlocksSingleStrand(String<uint64_t> & cords, String<UPair> & str_ends_p, 
    StringSet<String<UPair> > & cords_chains, ChainScoreMetric & chn_score, 
     uint64_t read_len, uint thd_init_cord_score)
{
    String<int> str_ends_p_score;
    resize(str_ends_p_score, length(str_ends_p));
    int strand = chn_score.chn_score_parms.chn_block_strand;
    if (strand)
    {
        std::sort (begin(str_ends_p), end(str_ends_p), [&cords, &read_len](UPair & a, UPair &b){
        uint64_t y1,y2;
        y1 = !get_cord_strand(cords[a.first]) ? read_len - 1 - get_cord_y(cords[a.second - 1]) :
            get_cord_y(cords[a.first]);
        y2 = !get_cord_strand(cords[b.first]) ? read_len - 1 - get_cord_y(cords[b.second - 1]) :
            get_cord_y(cords[b.first]);
        return y1 > y2;
        });
    }
    else
    {
        std::sort (begin(str_ends_p), end(str_ends_p), [&cords, &read_len](UPair & a, UPair &b){
        uint64_t y1,y2;
        y1 = get_cord_strand(cords[a.first]) ? read_len - 1 - get_cord_y(cords[a.second - 1]) :
            get_cord_y(cords[a.first]);
        y2 = get_cord_strand(cords[b.first]) ? read_len - 1 - get_cord_y(cords[b.second - 1]) :
            get_cord_y(cords[b.first]);
        return y1 > y2;
        });
    } 
    for (unsigned i = 0; i < length(str_ends_p_score); i++) 
    {
        //init the score as the length of the blocks
        str_ends_p_score[i] = (str_ends_p[i].second - str_ends_p[i].first) * thd_init_cord_score;
    }
    int thd_best_n1 = 3; //unlimited
    chainBlocksBase(cords_chains, cords, str_ends_p, str_ends_p_score, read_len, chn_score, thd_best_n1, 0);
}
/*
 * Return best strand: 0 @cords_chains1, 1:@cords_chains2
 */
int getChainBlocksBestStrand(StringSet<String<UPair> > & cords_chains1, 
                             StringSet<String<UPair> > & cords_chains2)
{
    String <int> lens1, lens2; //@lens[i] = length sum of @chains_cords[k], k<=i
    resize(lens1, length(cords_chains1));
    resize(lens2, length(cords_chains2));
    if (!empty(cords_chains1))
    {
        for (int i = 0; i < length(cords_chains1); i++)
        {
            lens1[i] = i == 0 ? 0 : lens2[i - 1];
            for (int j = 0; j < length(cords_chains1[i]); j++)
            {
                lens1[i] += cords_chains1[i][j].second - cords_chains1[i][j].first;
                dout << "gbs11" << i << cords_chains1[i][j].first << cords_chains1[i][j].second << "\n";
            }
        }
    }  
    if (!empty(cords_chains2))
    {
        for (int i = 0; i < length(cords_chains2); i++)
        {
            lens2[i] = i == 0 ? 0 : lens2[i - 1];
            for (int j = 0; j < length(cords_chains2[i]); j++)
            {
                lens2[i] += cords_chains2[i][j].second - cords_chains2[i][j].first;
                dout << "gbs12" << cords_chains2[i][j].first << cords_chains2[i][j].second << "\n";
            }
        }
    }
    for (int i = 0; i < std::min(length(lens1), length(lens2)); i++)
    {
        dout << "gbs" << length(cords_chains2) << lens1[i] << lens2[i] << "\n";
        if (lens1[i] < lens2[i])
        {
            return 1;
        }
        else if (lens1[i] > lens2[i])
        {
            return 0;
        }
    }
    return 0; 
}
/*
 * Revert order of blocks in @cords_chains if strands is different from @strand
 */
int revertChainBlockStrand(StringSet<String<UPair> > & cords_chains,
                           String<uint64_t> & cords,
                           int strand,
                           uint64_t read_len)
{
    uint64_t swap_str = 0;
    uint64_t f_strand = strand ? 1 : 0;
    UPair cordy_pair_pre;
    UPair cordy_pair;
    dout << "cbrs" << f_strand << "\n";
    for (unsigned i = 0; i < length(cords_chains); i++) 
    {
        appendValue(cords_chains[i], UPair(0,0));
        uint64_t strand_pre = 0;
        uint64_t strand_this = 0;
        for (unsigned j = 0; j < length(cords_chains[i]); j++) 
        {
            if (j == length(cords_chains[i]) - 1 || 
                get_cord_strand(cords[cords_chains[i][j].first]) == f_strand)
            {
                strand_this = 0;
            }
            else 
            {
                strand_this = 1;
            }
            if (strand_this && !strand_pre) 
            {
                swap_str = j;
            }
            if (!strand_this && strand_pre)
            {
                for (unsigned k = swap_str; k < (swap_str + j)/ 2; k++)
                {
                    std::swap(cords_chains[i][k], cords_chains[i][swap_str + j - 1 - k]);
                }
            }
            strand_pre = strand_this;
        }
        resize(cords_chains[i], length(cords_chains[i]) - 1);
    }
    return 0;
}
/*
 * chain blocks of @cords
   start and end @cords of each block is specified in @str_ends_p
 */
int chainBlocksCords(String<uint64_t> & cords, 
                     String<UPair> & str_ends_p, 
                     ChainScoreMetric & chn_score,  
                     uint64_t read_len,  
                     uint thd_init_cord_score, 
                     uint64_t thd_major_limit,
                     void (*unsetEndFunc)(uint64_t &), 
                     void(*setEndFunc)(uint64_t &), 
                     int f_header)
{
    StringSet<String<UPair> > cords_chains1;
    StringSet<String<UPair> > cords_chains2;
    String<UPair> str_ends_p1(str_ends_p);
    String<UPair> str_ends_p2(str_ends_p);
    chn_score.chn_score_parms.chn_block_strand = 0;
    dout << "cbs2" << length(str_ends_p1) << "\n";
        print_cords(cords, "cbs1");
    chainBlocksSingleStrand(cords, str_ends_p1, cords_chains1, chn_score, read_len, thd_init_cord_score);
    chn_score.chn_score_parms.chn_block_strand = 1;
    chainBlocksSingleStrand(cords, str_ends_p2, cords_chains2, chn_score, read_len, thd_init_cord_score);
    int best_strand = getChainBlocksBestStrand(cords_chains1, cords_chains2);
    if (best_strand == 0)
    {
        str_ends_p = str_ends_p1;
        revertChainBlockStrand(cords_chains1, cords, best_strand, read_len);
        _filterBlocksCords (cords_chains1, cords, read_len, thd_major_limit, unsetEndFunc, 
        setEndFunc, f_header);
        print_cords(cords, "cbs4");
    }
    else
    {
        str_ends_p = str_ends_p2;
        revertChainBlockStrand(cords_chains2, cords, best_strand, read_len);
        _filterBlocksCords (cords_chains2, cords, read_len, thd_major_limit, unsetEndFunc, 
        setEndFunc, f_header);  
    }
    print_cords(cords, "cbs3");
    return 0;
}

class NumericalScore
{
    float erf_num[32];
public:
    NumericalScore();
    float erf(float);
};
/*
 * Table of numerical approximation of error function
 */
NumericalScore::NumericalScore():
 erf_num{
//0,   0.02, 0.04, 0.06, 0.08, 0.1, 
//0.2, 0.3,  0.4,  0.5,  0.6,  0.7, 
//0.8, 0.9,  1,    1.1,  1.2,  1.3, 
//1.4, 1.5,  1.6,  1.7,  1.8,  1.9, 
//2.0, 2.1,  2.2,  2.3,  2.4,  2.5
//>2.5
0,           0.022564575, 0.045111106, 0.067621594, 0.090078126, 0.112462916,
0.222702589, 0.328626759, 0.428392355, 0.520499878, 0.603856091, 0.677801194,
0.742100965, 0.796908212, 0.842700793, 0.88020507,  0.910313978, 0.934007945,
0.95228512,  0.966105146, 0.976348383, 0.983790459, 0.989090502, 0.992790429,
0.995322265, 0.997020533, 0.998137154, 0.998856823, 0.999311486, 0.999593048,
1} 
{}
/*
 * Return numerical approximation of error function
 * val \in [-3.5, 3.5], otherwise return 1;
 */
float NumericalScore::erf(float val)
{
    float abs_val = val < 0 ? -val : val;
    float score = 0;
    if (abs_val > 2.5)
    {
        score = 1;  
    }
    else if (abs_val < 0.1) 
    {
        unsigned i = abs_val / 0.02;
        score = (erf_num[i] + erf_num[i + 1]) * 0.5;
    }
    else 
    {
        unsigned i =(5 + (abs_val - 0.1) / 0.1);
        score = (erf_num[i] + erf_num[i + 1]) * 0.5;
    }
        dout << "erf" << abs_val << score << "\n";
    return val < 0 ? -score:score;
}
NumericalScore num_score;
/*
 * return cdf function of normal function
 */
float cdfN(float val, float mean, float var)
{
    return (1 + num_score.erf((val - mean) / (var * 1.414))) * 0.5;
}
/*
 * Return probability of variants 
 * @strand: {0,1} same strand:=0; diff..:=1
 */
float variantsProb(int strand, int64_t dx, int64_t dy)
{
    int64_t da = dx - dy;
    float p = 1;  
    if (strand) //inv
    {
        p = 0.5;
    }
    if (da < -std::max(dx / 4, int64_t(50))) //ins, dup
    {
        if (dx > -50) //ins
        {
            p = 0.5; 
        }
        else 
        {
            p = 0.25;
        }
    }
    else if (da > std::max(dy / 4, int64_t(50))) //del
    {
        p = 0.5;
    }
    return p;
}
/*
 * Return score of chaining blocks of [cord11, cord12), and [cord21, cord22)
 */
int getChainBlocksScore1(uint64_t const & cord11, uint64_t const & cord12, 
    uint64_t const & cord21, uint64_t const & cord22, 
    uint64_t const & read_len, ChainScoreParms & chn_sc_parms)
{
    int64_t dx, dy;
    //int f_type = getForwardChainDxDy(cord11, cord12, cord21, cord22, read_len, dx, dy);
    int f_type = getChainBlockDxDy(cord11, cord12, cord21, cord22, read_len, chn_sc_parms.chn_block_strand, dx, dy);
    if (dy < -80)
    {
        return INT_MIN;
    }
    int64_t d = std::max(std::min(dx, dy), int64_t(0));
    float p_0 = 1 - cdfN(d, chn_sc_parms.mean_d, chn_sc_parms.var_d); //p_0 = probability of dist > d
    float p = variantsProb(f_type ? 1 : 0, dx, dy) * p_0;
    int score = p * 100; //100: float accuracy 0.01 
    dout << "gcbs1x" << 0 << cdfN(0.025, 0, 1) << cdfN(-0.025, 0, 1) << score << p << d << p_0 << "\n";

    return score;
}

//End all mapper module
