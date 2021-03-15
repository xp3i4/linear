#ifndef LINEAR_HEADER_ALIGN_BANDS_H
#define LINEAR_HEADER_ALIGN_BANDS_H
#include <seqan/align.h>
#include "f_io.h"
struct MergeCordsBandsParm 
{
    unsigned thd_search_depth;
    
    MergeCordsBandsParm();
};

int mergeCordsBands(String<uint64_t> & cords_str,
                    String<uint64_t> & cords_end,
                    String<uint64_t> & bands_lower,
                    String<uint64_t> & bands_upper,
                    MergeCordsBandsParm & parm);
#endif