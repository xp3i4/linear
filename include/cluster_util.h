#ifndef LINEAR_HEADER_CLUSTER_UTIL_H
#define LINEAR_HEADER_CLUSTER_UTIL_H
#include <seqan/sequence.h>
#include "base.h"
#include "cords.h"
using namespace seqan;
struct ChainScoreParms
{

};

int getApxChainScore(uint64_t const & anchor1, uint64_t const & anchor2, ChainScoreParms & chn_sc_parms);
int getApxChainScore2(uint64_t const & cord11, uint64_t const & cord12, uint64_t const & cord21, uint64_t const & cord22, uint64_t const & read_len, ChainScoreParms & chn_sc_parms);
int getApxChainScore3(uint64_t const & cord11, uint64_t const & cord12, uint64_t const & cord21, uint64_t const & cord22, uint64_t const & read_len, ChainScoreParms & chn_sc_parms);

//Chainning Score metric wrapper: including a score function with corresponding parms.
struct ChainScoreMetric 
{
    int thd_abort_score; //lower bound of average chain score per anchor
    ChainScoreParms chn_score_parms;

    int (*getScore) (uint64_t const &, uint64_t const &, ChainScoreParms &);
    int (*getScore2)(uint64_t const &, uint64_t const &, uint64_t const &, uint64_t const &, 
        uint64_t const & read_len, ChainScoreParms &);
    int getAbortScore();
    ChainScoreMetric();
    ChainScoreMetric(int abort_score, int (*scoreFunc) (uint64_t const &, uint64_t const &, ChainScoreParms &));
    ChainScoreMetric(int abort_score, int (*scoreFunc)(uint64_t const &, uint64_t const &, uint64_t const &, uint64_t const &, uint64_t const & read_len, ChainScoreParms &));
};

struct ChainsRecord
{
    int score; //chain score 
    int score2; //copy of score 
    int len;
    int p2anchor;
    int root_ptr;
    int f_leaf;

    int isLeaf();
};
int chainAnchorsHits(String<uint64_t> & anchors, String<uint64_t> & hits, String<int> & hits_chains_score);
int chainBlocksHits(String<uint64_t> & hits, String<UPair> & str_ends_p, String<int> & str_ends_p_score, uint64_t read_len);
int getBestChains(String<uint64_t> & anchor, String<ChainsRecord> & chains,
                  int (*getScore) (uint64_t const &, uint64_t const &));

int chainAnchorsBase(String<uint64_t> &, StringSet<String<uint64_t> > &, String<int> &, uint, uint, uint,
 uint64_t, int, ChainScoreMetric &, uint64_t (*get_anchor_x)(uint64_t));
int getForwarChainDxDy(uint64_t const & cord11, uint64_t const & cord12, uint64_t const & cord21, uint64_t const & cord22, uint64_t const & read_len, int64_t & dx, int64_t & dy);
int chainBlocksCords(String<uint64_t> & cords, String<UPair> & str_ends_p, ChainScoreMetric & chn_score,
   uint64_t read_len, uint thd_init_cord_score, uint64_t thd_major_limit, 
   void (*unsetEndFunc)(uint64_t &), void (*setEndFunc)(uint64_t &), int f_header);
#endif