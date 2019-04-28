#include <seqan/sequence.h>
#include "shape_extend.h" 
#include "pmpfinder.h"
#include "index_util.h"
using namespace seqan;

unsigned dshape_len = 14;
DIndex::DIndex():
    shape(dshape_len)
{}
DIndex::DIndex (unsigned len):
    shape(len)
{}
LShape & DIndex::getShape()
{
    return shape;
}
int DIndex::fullSize()
{
    return (1 << shape.weight << shape.weight) + 1;
}
String<int64_t> & DIndex::getDir()
{
    return dir;
}
String<int64_t> & DIndex::getHs()
{
    return hs;
}

int createDIndex(StringSet<String<Dna5> > & seqs, 
                 DIndex & index, 
                 int64_t thd_min_step, 
                 int64_t thd_max_step)
{
    double t = sysTime();
    LShape & shape = index.getShape();
    String<int64_t> tmp;
    String<int64_t> & dir = index.getDir();
    String<int64_t> & hs = index.getHs();
    resize (index.getDir(), index.fullSize(), 0);
    int64_t preVal = ~0;
    int64_t last_j = 0;
    for (int64_t i = 0; i < length(seqs); i++)
    {
        int64_t count = 0;
        hashInit (shape, begin(seqs[i]));
        for (int64_t j = 0; j < length(seqs[i]) - shape.span; j++)
        {
            hashNexth(shape, begin(seqs[i]) + j);
            if (++count > thd_min_step)
            {
                hashNextX(shape, begin(seqs[i]) + j);
                if (preVal != shape.XValue || j - last_j > thd_max_step)
                {
                    ++dir[shape.XValue];
                    preVal = shape.XValue;
                    last_j = j;
                }
                count = 0;
            }
        }
    }
    int64_t sum = 0;
    for (int64_t i = 0; i < length(dir); i++)
    {
        sum += dir[i];
        dir[i] = sum - dir[i];
    }
    last_j = 0;
    int64_t EmptyVal = create_cord(length(seqs),0,0,0); 
    //make sure genomeid >= length(seqs) and cord y be 0! y points to next empty.
    resize (hs, sum, EmptyVal);
    for (int64_t i = 0; i < length(seqs); i++)
    {
        int64_t count = 0; 
        hashInit (shape, begin(seqs[i]));
        for (int64_t j = 0; j < length(seqs[i]) - shape.span; j++)
        {
            hashNexth (shape, begin(seqs[i]) + j);
            if (++count > thd_min_step)
            {
                hashNextX (shape, begin(seqs[i]) + j);
                //std::cout << "cx4 " << j << " " << shape.XValue << " " << dir[shape.XValue] << "\n";
                if (preVal != shape.XValue || j - last_j > thd_max_step)
                {
                    int k = dir[shape.XValue];
                    k += get_cord_y (hs[k]);
                    hs[k] = create_cord(i, j, 0, shape.strand);
                    hs[dir[shape.XValue]] = _DefaultCord.shift(hs[dir[shape.XValue]], 0, 1); //y++;
                    //std::cout << shape.XValue << " " << k <<"xxx3\n";

                    //hs[k]++;
                    preVal = shape.XValue;
                    last_j = j;

                } 
                count = 0;
            }
        }  
    }
    std::cout << "createDIndex " << sysTime() - t << "\n";
}

int64_t queryHsStr(DIndex & index, int64_t xval)
{
    return index.getDir()[xval];
}
int64_t queryHsEnd(DIndex & index, int64_t xval)
{
    return index.getDir()[xval + 1];
}
