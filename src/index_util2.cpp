#include <seqan/sequence.h>
#include "shape_extend.h" 
#include "pmpfinder.h"
#include "index_util.h"
#include "base.h"
using namespace seqan;

int const typeDIx = 1;
int const typeHIx = 2;

unsigned dshape_len = 22;
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
String<int> & DIndex::getDir()
{
    return dir;
}
String<int64_t> & DIndex::getHs()
{
    return hs;
}
int createDIndex_serial(StringSet<String<Dna5> > & seqs, 
                        DIndex & index, 
                        int64_t thd_min_step, 
                        int64_t thd_max_step)
{
    double t = sysTime();
    LShape & shape = index.getShape();
    String<int> & dir = index.getDir();
    String<int64_t> & hs = index.getHs();
    resize (index.getDir(), index.fullSize(), 0);
    double t2 = sysTime();
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
                if (preVal != shape.XValue || j - last_j > thd_max_step)
                {
                    int k = dir[shape.XValue];
                    k += get_cord_y (hs[k]);
                    hs[k] = create_cord(i, j, 0, shape.strand);
                    hs[dir[shape.XValue]] = shift_cord(hs[dir[shape.XValue]], 0, 1); //y++;
                    //std::cout << shape.XValue << " " << k <<"xxx3\n";

                    //hs[k]++;
                    preVal = shape.XValue;
                    last_j = j;

                } 
                count = 0;
            }
        }  
    }
    std::cout << "createDIndex " << sysTime() - t << " " << sysTime() - t2 << "\n";
}

int createDIndex(StringSet<String<Dna5> > & seqs, 
                 DIndex & index, 
                 int64_t thd_min_step, 
                 int64_t thd_max_step,
                 int64_t thd_omit_block,
                 unsigned threads
                )
{
    serr.print_message(">>Index::initing ", 0, 2, std::cerr);
    double t = sysTime();
    LShape & t_shape = index.getShape();
    String<int> & dir = index.getDir();
    String<int64_t> & hs = index.getHs();
    resize (index.getDir(), index.fullSize(), 0);
    double t2 = sysTime();
    dout << threads << "threads\n"; 
    for (int64_t i = 0; i < length(seqs); i++)
    {
        String<int64_t> t_blocks;
        for (int j = 0; j < threads; j++)
        {
            appendValue(t_blocks, length(seqs[i]) / threads * j); 
        }
        appendValue (t_blocks, length(seqs[i]) - t_shape.span);
        dout << t_blocks << "cdx2\n";
    #pragma omp parallel
    {
        unsigned t_id = omp_get_thread_num();
        int64_t t_str = t_blocks[t_id];
        int64_t t_end = t_blocks[t_id + 1];
        int64_t preVal = ~0;
        int64_t last_j = t_str - 1;
        int64_t count = 0;
        LShape shape = t_shape;
        hashInit (shape, begin(seqs[i]) + t_str);
        dout << "cdx1 " << t_str<< t_end <<"\n";
        for (int64_t j = t_str; j < t_end; j++)
        {
            hashNexth(shape, begin(seqs[i]) + j);
            if (++count > thd_min_step)
            {
                hashNextX(shape, begin(seqs[i]) + j);
                if (preVal != shape.XValue || j - last_j > thd_max_step)
                {
                    atomicInc(dir[shape.XValue]);
                    //dout << get_cord_y(atomicInc(dir[shape.XValue])) << get_cord_y(dir[shape.XValue]) << "xxx1";
                    preVal = shape.XValue;
                    last_j = j;
                }
                count = 0;
            }
        }
    }
    }
    //double t4 = sysTime();
    int64_t sum = 0;
    for (int64_t i = 0; i < length(dir); i++)
    {
        if (dir[i] > thd_omit_block)
        {
            dir[i] = 0;
        }
        sum += dir[i];
        dir[i] = sum - dir[i];
    }
    //dout << sysTime() - t4 << "x5\n";
    int64_t EmptyVal = create_cord(length(seqs),0,0,0); 
    //make sure genomeid >= length(seqs) and cord y be 0! y points to next empty.
    resize (hs, sum, EmptyVal);
    serr.print_message("--Index::inite   ", 0, 1, std::cerr);
    serr.print_message(">>Index::hashing", 0, 2, std::cerr);
    for (int64_t i = 0; i < length(seqs); i++)
    {
        String<int64_t> t_blocks;
        for (int j = 0; j < threads; j++)
        {
            appendValue(t_blocks, length(seqs[i]) / threads * j); 
        }
        appendValue (t_blocks, length(seqs[i]) - t_shape.span);
        dout << "cx22"<<t_blocks << "\n";
    #pragma omp parallel
    {   
        unsigned t_id = omp_get_thread_num();
        int64_t t_str = t_blocks[t_id];
        int64_t t_end = t_blocks[t_id + 1];
        int64_t preVal = ~0;
        int64_t last_j = t_str - 1;
        int64_t count = 0;
        LShape shape = t_shape;
        hashInit (shape, begin(seqs[i]) + t_str);
        dout << "strd " << t_str << t_end << "\n";
        for (int64_t j = t_str; j < t_end; j++)
        {
            hashNexth (shape, begin(seqs[i]) + j);
            if (++count > thd_min_step)
            {
                hashNextX (shape, begin(seqs[i]) + j);
                if (preVal != shape.XValue || j - last_j > thd_max_step)
                {
                    if (dir[shape.XValue + 1] - dir[shape.XValue ])
                    {
                        int64_t k = dir[shape.XValue];
                        k += get_cord_y (atomic_inc_cord_y(hs[k])) - 1;
                        hs[k] = create_cord(i, j, 0, shape.strand);
                        hs[dir[shape.XValue]] = shift_cord(hs[dir[shape.XValue]], 0, 1); //y++;
                        preVal = shape.XValue;
                        last_j = j;
                    }
                } 
                count = 0;
            }
        }  
    }
    }
    std::cout << "createDIndex " << sysTime() - t << " " << sysTime() - t2 << "\n";
    serr.print_message("Index::hash        ", 2, 1, std::cerr);
    serr.print_message("End createing index", 2, 1, std::cerr);
}

int64_t queryHsStr(DIndex & index, int64_t xval)
{
    return index.getDir()[xval];
}
int64_t queryHsEnd(DIndex & index, int64_t xval)
{
    return index.getDir()[xval + 1];
}

int IndexDynamic::isHIndex()
{
    return typeIx == typeHIx;
}

int IndexDynamic::isDIndex()
{
    return typeIx == typeDIx;
}

void IndexDynamic::setHIndex()
{
    typeIx = typeHIx;
}

void IndexDynamic::setDIndex()
{
    typeIx = typeDIx;
}

IndexDynamic::IndexDynamic(StringSet<String<Dna5> > & seqs):hindex(seqs)
{}

bool createIndexDynamic(StringSet<String<Dna5> > & seqs, IndexDynamic & index, unsigned threads, bool efficient)
{
    int64_t thd_min_step = 4;
    int64_t thd_max_step = 10;
    int64_t thd_omit_block = 50;
    if (index.isDIndex())
    {
        std::cout << "cidx\n";
        //TODO::parm wrapping 
        return createDIndex(seqs, index.dindex, 
                            thd_min_step, 
                            thd_max_step, 
                            thd_omit_block,
                            threads);
//        return createDIndex_serial(seqs, index.dindex, 4, 10);
    }
    else if (index.isHIndex())
    {
        return createHIndex(seqs, index.hindex, threads, efficient);
    }

}

