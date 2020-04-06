#include <iostream>
#include "base.h"
#include "cords.h"
#include "shape_extend.h"
#include "index_util.h"
#include "cluster_util.h"
#include "pmpfinder.h"

using namespace seqan;
using std::cout;
using std::endl;

int const typeFeatures1_16 = 0;
int const typeFeatures1_32 = 1;
int const typeFeatures2_48 = 2;

int FeaturesDynamic::isFs1_16()
{
    return fs_type == typeFeatures1_16;
}
int FeaturesDynamic::isFs1_32()
{
    return fs_type == typeFeatures1_32;
}
int FeaturesDynamic::isFs2_48()
{
    return fs_type == typeFeatures2_48;
}
void FeaturesDynamic::setFs1_16()
{
    fs_type = typeFeatures1_16;
}
void FeaturesDynamic::setFs1_32()
{
    fs_type = typeFeatures1_32;
}
void FeaturesDynamic::setFs2_48()
{
    fs_type = typeFeatures2_48;
}
void FeaturesDynamic::setFeatureType(int type)
{
    if (type == typeFeatures1_16)
    {
        setFs1_16();
    }
    else if (type == typeFeatures1_32)
    {
        setFs1_32();
    }
    else //default 
    {
        setFs2_48();
    }
}
ApxMapParm1_16 _apx_parm1_16;
ApxMapParm1_32 _apx_parm1_32;
ApxMapParm2_48 _apx_parm2_48;
int FeaturesDynamic::init(int type)
{
    setFeatureType(type);
    if (isFs1_16())
    {
        apx_parm1_16 = new ApxMapParm1_16();
    }
    if (isFs1_32())
    {
        apx_parm1_32 = new ApxMapParm1_32();
    }
    else if (isFs2_48())
    {
        apx_parm2_48 = new ApxMapParm2_48(); 
    }  
    return fs_type;
}
FeaturesDynamic::FeaturesDynamic(){}
FeaturesDynamic::FeaturesDynamic(int type)
{
    init(type);
}

void decInt96 (int96 & val, int96 & dval) // val -= dval
{
    val[0] -= dval[0];
    val[1] -= dval[1];
    val[2] -= dval[2];
}
void setInt96 (int96 & val, int96 & dval) // val = dval
{
    val[0] = dval[0];
    val[1] = dval[1];
    val[2] = dval[2];
}
void incInt96 (int96 & val, int96 & dval) // val += dval
{
    val[0] += dval[0];
    val[1] += dval[1];
    val[2] += dval[2];
}
void printInt96(int96 val, CharString header)
{
    std::cout << header << " ";
    int sum = 0;
    for (int i = 0; i < 3; i++)
    {
        int v = val[i];
        for (int ii = 0; ii < 5; ii++)
        {
            sum += v & 63;
            std::cout << (v & 63) << " ";
            v >>= 6;
        }
        std::cout << "| ";
    }
    std::cout << sum << "\n";
}

//Map "N" to a larger number > 64, such that (1LL << ordN_(N)) = 0.
inline int64_t ordN_(Dna5 & a)
{
    int64_t res = ordValue(a);
    return (res == 4) ? 100LL : res;
}

/*===================================================
=            Approximate mapping section            =
===================================================*/
//ApxMapParm
//variable marked by * is allowed to be modified 
//variable of $ must be changed correspondingly when * is modified
//varibel of const is not allowed to change.
unsigned window_size;
ApxMapParmBase::ApxMapParmBase (float v1,
                                unsigned v2,
                                unsigned v3,
                                unsigned v4, 
                                unsigned v42,
                                unsigned v5,
                                unsigned v6,
                                unsigned v7,
                                unsigned v8) :
    band_width(v1),
    cell_size(v2),
    cell_num(v3),
    windowThreshold(v4),
    windowThresholdReject(v42),
    windowSize(cell_size * cell_num), 
    windowDelta(windowSize * (1 - 2 * band_width)),
    sup(cell_num),
    med(ceil((1 - band_width) * cell_num)),
    inf(ceil((1 - 2 * band_width) * cell_num)),
    scpt_step(v5), //const
    scpt_bit(v6),  //const 2 ^ scpt_bit = scpt_step
    scpt_size(v7),
    scpt_num(windowSize / scpt_size), //const
    scpt_int_step (scpt_size / scpt_step), //const
    abort_score(v8)
{
    window_size = windowSize;
}

ApxMapParm1_16::ApxMapParm1_16():
    ApxMapParmBase(0.25, 16, 12, 60, 80, 16, 4, 16, 1000), //todo::reject score 50 needs test
    //----scprit parm----
    scpt_len(5), //const 2 ^ scpt_len == 32
    scpt_len2(scpt_len << 1), //const
    scriptMask((1 << scpt_len) - 1), //const
    scriptMask2(scriptMask << scpt_len), //const
    scptCount{1, 1 << scpt_len, 1 << (scpt_len * 2), 0, 0} //const
{
    window_size = windowSize;
}

ApxMapParm1_32::ApxMapParm1_32():
    ApxMapParmBase(0.25, 16, 12, 36, 50, 16, 4, 32, 1000), //todo::reject score 50 needs test
    //----scprit parm----
    scpt_len(5), //const 
    scpt_len2(scpt_len << 1), //const
    scriptMask((1 << scpt_len) - 1), //const
    scriptMask2(scriptMask << scpt_len), //const
    scptCount{1, 1 << scpt_len, 1 << (scpt_len * 2), 0, 0} //const
{
    window_size = windowSize;
}

ApxMapParm2_48::ApxMapParm2_48():
    ApxMapParmBase(0.25, 16, 6, 36, 50, 16, 4, 48, 1000)
{
    window_size = windowSize;
}
unsigned getFeatureWindowDelta(FeaturesDynamic & fs)
{
    if (fs.isFs2_48())
    {
        return fs.apx_parm2_48->windowDelta;
    }
    else if (fs.isFs1_32())
    {
        return fs.apx_parm1_32->windowDelta;
    }
    else if (fs.isFs1_16())
    {
        return fs.apx_parm1_16->windowDelta;
    }
}
unsigned getFeatureWindowDelta(StringSet<FeaturesDynamic> & fss)
{
    if (!empty(fss))
    {
        return getFeatureWindowDelta(fss[0]);
    }
    else
    {
        //return _apx_parm_base.windowDelta;
        return 0;
    }
}
unsigned getFeatureWindowSize(FeaturesDynamic & fs)
{
    if (fs.isFs2_48())
    {
        return fs.apx_parm2_48->windowSize;
    }
    else if (fs.isFs1_32())
    {
        return fs.apx_parm1_32->windowSize;
    }
    else if (fs.isFs1_16())
    {
        return fs.apx_parm1_16->windowSize;
    }
}
unsigned getFeatureWindowSize(StringSet<FeaturesDynamic> & fss)
{
    if (!empty(fss))
    {
        return getFeatureWindowSize(fss[0]);
    }
    else
    {
        //return _apx_parm_base.windowSize;
        return 0;
    }
}
unsigned getWindowThreshold(FeaturesDynamic & fs)
{
    if (fs.isFs2_48())
    {
        return fs.apx_parm2_48->windowThreshold;
    }
    else if (fs.isFs1_32())
    {
        return fs.apx_parm1_32->windowThreshold;
    }
    else if (fs.isFs1_16())
    {
        return fs.apx_parm1_16->windowThreshold;
    }
    //return _apx_parm_base.windowThreshold;
    return 1000;
}
unsigned getWindowThreshold(StringSet<FeaturesDynamic> & fss)
{
    if (!empty(fss))
    {
        return getWindowThreshold(fss[0]);
    }
    else
    {
        return 1000;
    }
}

unsigned getWindowThresholdReject(FeaturesDynamic & fs)
{
    if (fs.isFs2_48())
    {
        return fs.apx_parm2_48->windowThresholdReject;
    }
    else if (fs.isFs1_32())
    {
        return fs.apx_parm1_32->windowThresholdReject;
    }
    else if (fs.isFs1_16())
    {
        return fs.apx_parm1_16->windowThresholdReject;
    }
    //return _apx_parm_base.windowThreshold;
    return 1000; //reject all
}
unsigned getWindowThresholdReject(StringSet<FeaturesDynamic> & fss)
{
    if (!empty(fss))
    {
        return getWindowThresholdReject(fss[0]);
    }
    else
    {
        return 0; //reject all
    }
}


/*----------  Script encoding type I 1_32  ----------*/
//script:16 bits short int, contains 3 segments
unsigned __scriptDist16_3(short const & s1, short const & s2, 
                         short const & mask, unsigned const & len1, unsigned const & len2)
{
    return  std::abs((s1 & mask) - (s2 & mask)) +
            std::abs(((s1 >> len1) & mask) - ((s2 >> len1) & mask)) + 
            std::abs((s1 >> len2) - (s2 >> len2));
}
unsigned _scriptDist1_32(short const & s1, short const & s2, ApxMapParm1_32 & parm)
{
    return __scriptDist16_3(s1, s2, parm.scriptMask, parm.scpt_len, parm.scpt_len2);
}
unsigned _windowDist1_32(Iterator<String<short> >::Type const & it1, 
                         Iterator<String<short> >::Type const & it2,
                         ApxMapParm1_32 & parm)
{
    unsigned dist = 0;
    for (unsigned i = 0; i < parm.scpt_num * parm.scpt_int_step; i += parm.scpt_int_step) 
    {
        dist += _scriptDist1_32(*(it1 + i), *(it2 + i), parm);
    }
    return dist;
}
void createFeatures1_32(TIter5 const & itBegin, TIter5 const & itEnd, String<short> & f, ApxMapParm1_32 & parm)
{
    unsigned next = 1;
    unsigned window = parm.scpt_size;
    resize (f, ((itEnd - itBegin - window) >> parm.scpt_bit) + 1);
    f[0] = 0;
    for (unsigned k = 0; k < window; k++)
    {
        f[0] += parm.scptCount[ordValue(*(itBegin + k))];
    }
    for (unsigned k = parm.scpt_step; k < itEnd - itBegin - window ; k += parm.scpt_step) 
    {
        f[next] = f[next - 1];
        for (unsigned j = k - parm.scpt_step; j < k; j++)
        {
            f[next] += parm.scptCount[ordValue(*(itBegin + j + window))] - parm.scptCount[ordValue(*(itBegin + j))];
        }
        next++;
    }
}
uint64_t parallelParm_Static(uint64_t range, 
                             unsigned threads, 
                             unsigned & thd_id, 
                             uint64_t & thd_begin, 
                             uint64_t & thd_end)
{
    uint64_t ChunkSize = range / threads;
    unsigned id = range - ChunkSize * threads;
    if (thd_id < id)
    {
        thd_begin = ++ChunkSize * thd_id; 
    }
    else
    {
       thd_begin = id + thd_id * ChunkSize;
    }
    thd_end = thd_begin + ChunkSize; 
    return ChunkSize;
}
void createFeatures1_32(TIter5 const & itBegin, TIter5 const & itEnd, String<short> & f, unsigned threads, ApxMapParm1_32 & parm)
{
    unsigned window = parm.scpt_size;
    resize (f, ((itEnd - itBegin - window) >> parm.scpt_bit) + 1);
#pragma omp parallel
{
    unsigned thd_id = omp_get_thread_num();
    uint64_t thd_begin;
    uint64_t thd_end;
    uint64_t range = (itEnd - itBegin - window - parm.scpt_step) / parm.scpt_step;
    parallelParm_Static(range, threads, 
                        thd_id,  thd_begin, thd_end);
    uint64_t next = thd_begin;
    thd_begin *= parm.scpt_step;
    thd_end *= parm.scpt_step;
    f[next] = 0;
    for (unsigned k = thd_begin; k < thd_begin + window; k++)
    {
        f[next] += parm.scptCount[ordValue(*(itBegin + k))];
    }
    next++;
    for (unsigned k = thd_begin + parm.scpt_step; k < thd_end ; k += parm.scpt_step) 
    {
        f[next] = f[next - 1];
        for (unsigned j = k - parm.scpt_step; j < k; j++)
        {
            f[next] += parm.scptCount[ordValue(*(itBegin + j + window))] - parm.scptCount[ordValue(*(itBegin + j))];
        }
        next++;
    }
}
}

/*----------  Script encoding type I.2 1_16 ----------*/
unsigned _scriptDist1_16(int const & s1, int const & s2, ApxMapParm1_16 & parm)
{
    return __scriptDist16_3(s1, s2, parm.scriptMask, parm.scpt_len, parm.scpt_len2);
}
unsigned _windowDist1_16(Iterator<String<short> >::Type const & it1, 
                         Iterator<String<short> >::Type const & it2,
                         ApxMapParm1_16 & parm)
{
    unsigned sum = 0;
    for (unsigned i = 0; i < parm.scpt_num * parm.scpt_int_step; i += parm.scpt_int_step) 
    {
        sum += _scriptDist1_16(*(it1 + i), *(it2 + i), parm);
    }
    return sum;
}
void createFeatures1_16(TIter5 const & itBegin, TIter5 const & itEnd, String<short> & f, ApxMapParm1_16 & parm)
{
    unsigned window = 16;
    resize (f, ((itEnd - itBegin - window) >> parm.scpt_bit) + 1);
    uint64_t next = 0;
    for (unsigned k = 0; k < itEnd - itBegin - window ; k += parm.scpt_step) 
    {
        f[next] = 0;
        for (unsigned j = k; j < k + window; j++)
        {
            f[next] += parm.scptCount[ordValue(*(itBegin + j))];
        }
        next++;
    }
}
void createFeatures1_16(TIter5 const & itBegin, TIter5 const & itEnd, String<short> & f, unsigned threads, ApxMapParm1_16 & parm)
{
    unsigned window = 16;
    resize (f, ((itEnd - itBegin - window) >> parm.scpt_bit) + 1);
#pragma omp parallel
{
    unsigned thd_id = omp_get_thread_num();
    uint64_t thd_begin;
    uint64_t thd_end;
    uint64_t range = (itEnd - itBegin - window - parm.scpt_step) / parm.scpt_step;
    parallelParm_Static(range, threads, 
                        thd_id,  thd_begin, thd_end);
    uint64_t next = thd_begin;
    thd_begin *= parm.scpt_step;
    thd_end *= parm.scpt_step;

    for (unsigned k = thd_begin; k < thd_end - window; k += parm.scpt_step) 
    {
        f[next] = 0;
        for (unsigned j = k; j < k + window; j++)
        {
            f[next] += parm.scptCount[ordValue(*(itBegin + j))];
        }
        next++;
    }
}
}

/*----------  Script encoding type II  ----------*/
/* TG TC TA GT GG int1  \
   GC GA CT CG CC int2  -- int96 for 48 bases script; Each 2mer has 6 bits
   CA AT AG AC AA int3  /
   features[i] for 48 bases. 
 */
//TODO!!! for 2mer > 32 occurrences out of boundary

int const window48 = 48;
int const max31 = 31;
int const mxu31 = (max31 << 24) + (max31 << 18) + (max31 << 12) +
                  (max31 << 6) + max31;
int64_t __scriptDist63_31(int const & s1, int const & s2)
{
    int d = s1 + mxu31 - s2;
    int mask = 63; 
    return std::abs((d >> 24 & mask) - max31) +
           std::abs((d >> 18 & mask) - max31) +
           std::abs((d >> 12 & mask) - max31) +
           std::abs((d >> 6 & mask) - max31) +
           std::abs((d & mask) - max31);
}
int64_t _scriptDist63_31(int const & s11, int const & s12, int const & s13,
                         int const & s21, int const & s22, int const & s23)
{
    return __scriptDist63_31(s11, s21) +
           __scriptDist63_31(s12, s22) +
           __scriptDist63_31(s13, s23);
}
//wrapper into 96 bit integer.
//Distance of one pair of scpirts(48mer of each)
int64_t _scriptDist63_31(int96 & s1, int96 & s2)
{
    return  _scriptDist63_31(s1[0], s1[1], s1[2],
                             s2[0], s2[1], s2[2]);
}
int64_t _windowDist2_48(Iterator<String<int96> >::Type it1,  //feature string iterator
                        Iterator<String<int96> >::Type it2,
                        ApxMapParm2_48 & parm)
{
    int64_t sum = 0;
    for (unsigned i = 0; i < parm.scpt_num * parm.scpt_int_step; i += parm.scpt_int_step) 
    {
        sum += _scriptDist63_31(*(it1 + i), *(it2 + i));
    }
    return sum;
}
/**
 * 1.units[n] = (i << 8 + k) maps n to bits of int96 (ith int, kth bit);
 * 2.NOTE::In gcc, x << y == x << (y|31) due to optimization.
 *   Namely if y > 31, do nothing. 
 *   Hence the infiN is set to 31 for TT, N* and *N 
 *   They are add to the 31bits and cut off the 30bits by & infi_mask30
 */
short infiN = 31; //for base '*N' 'N*' and 'TT' 
unsigned infi_mask30 = (1 << 31) - 1;
short units[25] = {
/*A**/ 0,             6,             12,            18,            infiN, 
/*C**/ 24,            (1 << 8) + 0,  (1 << 8) + 6,  (1 << 8) + 12, infiN,
/*G**/ (1 << 8) + 18, (1 << 8) + 24, (2 << 8) + 0,  (2 << 8) + 6,  infiN,
/*T**/ (2 << 8) + 12, (2 << 8) + 18, (2 << 8) +24,  infiN,         infiN,  
/*N**/ infiN,         infiN,         infiN,         infiN,         infiN};
void add2merInt96 (int96 & val, TIter5 it)
{
    unsigned ordV = ordValue(*it) * 5 + ordValue(*(it + 1)); // 5*ordV_i + ordV_{i+1}
    unsigned i = units[ordV] >> 8;
    unsigned addVal = (1 << (units[ordV] & 255)) & infi_mask30;
    val[i] += addVal;
}
int createFeatures2_48(TIter5 it_str, TIter5 it_end, String<int96> & f, ApxMapParm2_48 & parm)
{
    double t = sysTime();
    int addMod3[3] = {1, 2, 0};   // adMod3[i] = ++ i % 3;
    int96 zero96 = {0, 0, 0};
    std::vector<int96> buffer(3, zero96); //buffer of 3 cells in one script
    resize (f, (it_end - it_str - window48) / parm.scpt_step + 1); 
    setInt96(f[0], zero96);
    for (int i = 0; i < 3; i++) //init f[0]
    {
        for (int j = i << parm.scpt_bit; j < (i << parm.scpt_bit) + parm.scpt_step; j++)
        {
            add2merInt96(buffer[i], it_str + j);
        }
        incInt96(f[0], buffer[i]);
    }
    //std::cerr << " cf48 time " << sysTime() - t << "\n";
    int next = 1; //stream f[next]
    int ii = 0;
    for (int i = parm.scpt_step; i < it_end - it_str - window48 - 1; i += parm.scpt_step) 
    {
        setInt96(f[next], f[next - 1]);
        decInt96(f[next], buffer[ii]);
        setInt96(buffer[ii], zero96); //clear to 0;
        for (int j = i - parm.scpt_step + window48; j < i + window48; j++)
        {
            add2merInt96(buffer[ii], it_str + j);
        }
        incInt96(f[next], buffer[ii]);
        ii = addMod3[ii];
        next++;
    }
    return 0;
}
int createFeatures2_48(TIter5 it_str, TIter5 it_end, String<int96> & f, unsigned threads, ApxMapParm2_48 & parm)
{
    int window = window48;
    if (it_end - it_str < window)
    {
        return 0;
    }
    resize (f, ((it_end - it_str - window) >> parm.scpt_bit) + 1);
    int64_t range = (it_end - it_str - window) / parm.scpt_step + 1; //numer of windows 
    if (range < threads)
    {
        createFeatures2_48(it_str, it_end, f, parm);
        return 0;
    }
#pragma omp parallel
{
    unsigned thd_id = omp_get_thread_num();
    int64_t chunk_size = range / threads; 
    int64_t thd_begin = thd_id * (chunk_size + 1);
    unsigned id1 = range - chunk_size * threads;
    if (thd_id >= id1)
    {
        thd_begin = id1 + chunk_size * thd_id;
    }
    else
    {
        ++chunk_size;
    }
    int64_t thd_end = thd_begin + chunk_size;
    int64_t next = thd_begin;
    thd_begin *= parm.scpt_step;
    thd_end *= parm.scpt_step;

    int addMod3[3] = {1, 2, 0};   // adMod3[i] = ++ i % 3;
    int96 zero96 = {0, 0, 0};
    std::vector<int96> buffer(3, zero96); //buffer of 3 cells in one script
    setInt96(f[next], zero96);
    for (int i = 0; i < 3; i++) //init buffer and f[next]
    {
        int tmp = thd_begin + (i << parm.scpt_bit); 
        for (int j = tmp; j < tmp + parm.scpt_step; j++)
        {
            add2merInt96(buffer[i], it_str + j);
        }
        incInt96(f[next], buffer[i]);
    }
    int ii = 0;
    next++; //stream f[next]
    for (int i = thd_begin + parm.scpt_step; i < thd_end; i += parm.scpt_step) 
    {
        setInt96(f[next], f[next - 1]);
        decInt96(f[next], buffer[ii]);
        setInt96(buffer[ii], zero96); //clear to 0;
        for (int j = i - parm.scpt_step + window; j < i + window; j++)
        {
            add2merInt96(buffer[ii], it_str + j);
        }
        incInt96(f[next], buffer[ii]);
        ii = addMod3[ii];
        next++;
    }
}
return 0;
}

/*----------  Script encoding wrapper  ----------*/
unsigned __windowDist(FeaturesDynamic & f1,
                      FeaturesDynamic & f2,
                      uint64_t x1, uint64_t x2)
{
    
    if (f1.isFs2_48())
    {
        return _windowDist2_48 (begin(f1.fs2_48) + x1, begin(f2.fs2_48) + x2, *(f2.apx_parm2_48));
    }
    else if (f1.isFs1_16())
    {
        return _windowDist1_16 (begin(f1.fs1_16) + x1, begin(f2.fs1_16) + x2, *(f1.apx_parm1_16));
    }
    else if (f1.isFs1_32())
    {
        return _windowDist1_32 (begin(f1.fs1_32) + x1, begin(f2.fs1_32) + x2, *(f1.apx_parm1_32));
    }
}

//The wrapper is(only) used in the gap.cpp
//Do not call this function frequently since the condition branch will drain the performance.
//NOTE::boundary of features is checked in this function
unsigned _windowDist(FeaturesDynamic & f1,
                     FeaturesDynamic & f2,
                     uint64_t x1, uint64_t x2)
{
    if (f1.isFs2_48())
    {
        if (x1 < length(f1.fs2_48) && x2 < length(f2.fs2_48))
        {
            return _windowDist2_48 (begin(f1.fs2_48) + x1, begin(f2.fs2_48) + x2, *(f2.apx_parm2_48));
        }
        else
        {
            return f1.apx_parm2_48->abort_score;
        }
    }
    else if (f1.isFs1_16())
    {
        if (x1 < length(f1.fs1_16) && x2 < length(f2.fs1_16))
        {
            return _windowDist1_16 (begin(f1.fs1_16) + x1, begin(f2.fs1_16) + x2, * (f1.apx_parm1_16));
        }
        else
        {
            return f1.apx_parm1_16->abort_score;
        }
    }
    else if (f1.isFs1_32())
    {
        if (x1 < length(f1.fs1_32) && x2 < length(f2.fs1_32))
        {
            return _windowDist1_32 (begin(f1.fs1_32) + x1, begin(f2.fs1_32) + x2, *(f1.apx_parm1_32));
        }
        else
        {
            return f1.apx_parm1_32->abort_score;
        }
    }
}

int createFeatures(TIter5 it_str, TIter5 it_end, FeaturesDynamic & f)
{

    if (f.isFs2_48())
    {
        createFeatures2_48(it_str, it_end, f.fs2_48, *(f.apx_parm2_48));
    }
    else if (f.isFs1_16())
    {
        createFeatures1_16(it_str, it_end, f.fs1_16, *(f.apx_parm1_16));
    }
    else if (f.isFs1_32())
    {
        createFeatures1_32(it_str, it_end, f.fs1_32, *(f.apx_parm1_32));
    }
    return 0;
}

int createFeatures(TIter5 it_str, TIter5 it_end, FeaturesDynamic & f, unsigned threads)
{

    if (f.isFs2_48())
    {
        createFeatures2_48(it_str, it_end, f.fs2_48, threads, *(f.apx_parm2_48));
    }
    else if (f.isFs1_16())
    {
        createFeatures1_16(it_str, it_end, f.fs1_16, threads, *(f.apx_parm1_16));
    }
    else if (f.isFs1_32())
    {
        createFeatures1_32(it_str, it_end, f.fs1_32, threads, *(f.apx_parm1_32));
    }
    return 0;
}

//serial
int createFeatures(StringSet<String<Dna5> > & seq, 
                   StringSet<FeaturesDynamic> & f,
                   int feature_type)
{
    resize(f, length(seq));
    for (unsigned k = 0; k < length(seq); k++)
    {
        f[k].init(feature_type);
        createFeatures(begin(seq[k]), end(seq[k]), f[k]);
    }
    return 0;
}

//parallel
int createFeatures(StringSet<String<Dna5> > & seq, 
                   StringSet<FeaturesDynamic> & f, 
                   int feature_type,
                   unsigned threads)
{
    resize(f, length(seq));
    for (unsigned k = 0; k < length(seq); k++)
    {
        f[k].init(feature_type);
        createFeatures(begin(seq[k]), end(seq[k]), f[k], threads);
    }
    return 0;
}

/*----------  Dynamic programming of extending path (tiles)  ----------*/

//Template function is not used to speed up the compilation 
uint64_t previousWindow1_32(String<short> & f1, 
                            String<short> & f2, 
                            uint64_t cord,
                            float & score,
                            ApxMapParm1_32 & parm)
{
    uint64_t genomeId = get_cord_id(cord);
    uint64_t strand = get_cord_strand(cord);
    uint64_t x_suf = _DefaultCord.cord2Cell(get_cord_x(cord));
    uint64_t y_suf = _DefaultCord.cord2Cell(get_cord_y(cord));
    uint64_t x_min = 0;
    uint64_t y;
    uint64_t new_cord = 0;
    
    if (y_suf < parm.med || x_suf < parm.sup)
        return 0;
    else 
        y = y_suf - parm.med;

    unsigned min = ~0;
    for (uint64_t x = x_suf - parm.sup; x < x_suf - parm.inf; x += 1) 
    {
        unsigned tmp = _windowDist1_32(begin(f1) + y, begin(f2) + x, parm);
        if (tmp < min)
        {
            min = tmp;
            x_min = x;
        }
    }
    if (min > parm.windowThreshold)
        return 0;    
    else 
    {
        if ( x_suf - x_min > parm.med)
        {
            new_cord = _DefaultCord.createCord(create_id_x(genomeId, _DefaultCord.cell2Cord(x_suf - parm.med)),  _DefaultCord.cell2Cord(x_suf - x_min - parm.med + y), strand);
        }
        else
        {
            new_cord = _DefaultCord.createCord(create_id_x(genomeId, _DefaultCord.cell2Cord(x_min)), _DefaultCord.cell2Cord(y), strand);
        } 
    }
    score = min;
    return new_cord;
}

uint64_t previousWindow2_48(String<int96> & f1, 
                            String<int96> & f2, 
                            uint64_t cord,
                            float & score,
                            ApxMapParm2_48 & parm)
{
    uint64_t genomeId = get_cord_id(cord);
    uint64_t strand = get_cord_strand(cord);
    uint64_t x_suf = _DefaultCord.cord2Cell(get_cord_x(cord));
    uint64_t y_suf = _DefaultCord.cord2Cell(get_cord_y(cord));
    uint64_t x_min = 0;
    uint64_t y;
    uint64_t new_cord = 0;
    
    if (y_suf < parm.med || x_suf < parm.sup)
        return 0;
    else 
        y = y_suf - parm.med;

    unsigned min = ~0;
    for (uint64_t x = x_suf - parm.sup; x < x_suf - parm.inf; x += 1) 
    {
        unsigned tmp = _windowDist2_48(begin(f1) + y, begin(f2) + x, parm);
        if (tmp < min)
        {
            min = tmp;
            x_min = x;
        }
    }
    if (min > parm.windowThreshold)
        return 0;    
    else 
    {
        if ( x_suf - x_min > parm.med)
        {
            new_cord = _DefaultCord.createCord(create_id_x(genomeId, _DefaultCord.cell2Cord(x_suf - parm.med)),  _DefaultCord.cell2Cord(x_suf - x_min - parm.med + y), strand);
        }
        else
        {
            new_cord = _DefaultCord.createCord(create_id_x(genomeId, _DefaultCord.cell2Cord(x_min)), _DefaultCord.cell2Cord(y), strand);
        } 
    }
    score = min;
    return new_cord;
}
uint64_t previousWindow(FeaturesDynamic & f1, //read  
                        FeaturesDynamic & f2, 
                        uint64_t cord,
                        float & score)
{
    uint64_t genomeId = get_cord_id(cord);
    uint64_t strand = get_cord_strand(cord);
    uint64_t x_suf = _DefaultCord.cord2Cell(get_cord_x(cord));
    uint64_t y_suf = _DefaultCord.cord2Cell(get_cord_y(cord));
    uint64_t x_min = 0;
    uint64_t y;
    uint64_t new_cord = 0;
    ApxMapParmBase * parm;
    unsigned len1, len2; 

    if (f1.isFs2_48())
    {
        parm = f1.apx_parm2_48;
        len1 = length(f1.fs2_48);
        len2 = length(f2.fs2_48);
    }
    else if (f1.isFs1_16())
    {
        parm = f1.apx_parm1_16;
        len1 = length(f1.fs1_16);
        len2 = length(f2.fs1_16);
    }
    else if (f1.isFs1_32())
    {
        parm = f2.apx_parm1_32;
        len1 = length(f1.fs1_32);
        len2 = length(f2.fs1_32);
    }
    
    if (y_suf < parm->med || x_suf < parm->sup)
        return 0;
    else 
        y = y_suf - parm->med;
    unsigned min = ~0;
    for (uint64_t x = x_suf - parm->sup; x < x_suf - parm->inf; x += 1) 
    {
        unsigned tmp = __windowDist(f1, f2, y, x);
        if (tmp < min)
        {
            min = tmp;
            x_min = x;
        }
    }
    if (min > parm->windowThreshold)
        return 0;    
    else 
    {
        if (x_suf - x_min > parm->med)
        {
            new_cord = _DefaultCord.createCord(create_id_x(genomeId, _DefaultCord.cell2Cord(x_suf - parm->med)),  _DefaultCord.cell2Cord(x_suf - x_min - parm->med + y), strand);
        }
        else
        {
            new_cord = _DefaultCord.createCord(create_id_x(genomeId, _DefaultCord.cell2Cord(x_min)), _DefaultCord.cell2Cord(y), strand);
        } 
    }
    score = min;
    return new_cord;
}

uint64_t nextWindow1_32(String<short> & f1, 
                        String<short> & f2, 
                        uint64_t cord,
                        float & score,
                        ApxMapParm1_32 & parm)
{
    uint64_t genomeId = get_cord_id(cord);
    uint64_t strand = get_cord_strand(cord);
    uint64_t x_pre = _DefaultCord.cord2Cell(get_cord_x(cord));
    uint64_t y_pre = _DefaultCord.cord2Cell(get_cord_y(cord));
    uint64_t x_min = 0;
    uint64_t y;
    uint64_t new_cord = 0;
    unsigned min = ~0;
    
    if (y_pre + parm.sup * 2 > length(f1) || x_pre + parm.sup * 2> length(f2))
        return 0;
    else 
        y = y_pre + parm.med;

    for (uint64_t x = x_pre + parm.inf; x < x_pre + parm.sup; x += 1) 
    {
        unsigned tmp = _windowDist1_32(begin(f1) + y, begin(f2) + x, parm);
        if (tmp < min)
        {
            min = tmp;
            x_min = x;
        }
    }
    if (min > parm.windowThreshold)
    {
       return 0;
    }
    else 
    {
        if ( x_min - x_pre > parm.med)
        {
            new_cord = _DefaultCord.createCord(create_id_x(genomeId, _DefaultCord.cell2Cord(x_pre + parm.med)),  _DefaultCord.cell2Cord(x_pre + parm.med - x_min + y), strand);
        }
        else
        {
            new_cord = _DefaultCord.createCord(create_id_x(genomeId, _DefaultCord.cell2Cord(x_min)), _DefaultCord.cell2Cord(y), strand);
        }
    }
    score = min;
    return new_cord;
}

uint64_t nextWindow2_48(String<int96> & f1, //read  
                        String<int96> & f2, 
                        uint64_t cord,
                        float & score,
                        ApxMapParm2_48 & parm)
{
    uint64_t genomeId = get_cord_id(cord);
    uint64_t strand = get_cord_strand(cord);
    uint64_t x_pre = _DefaultCord.cord2Cell(get_cord_x(cord));
    uint64_t y_pre = _DefaultCord.cord2Cell(get_cord_y(cord));
    uint64_t x_min = 0;
    uint64_t y;
    uint64_t new_cord = 0;
    unsigned min = ~0;
    
    if (y_pre + parm.sup * 2 > length(f1) || x_pre + parm.sup * 2 > length(f2))
        return 0;
    else 
        y = y_pre + parm.med;
    
    for (uint64_t x = x_pre + parm.inf; x < x_pre + parm.sup; x += 1) 
    {
        unsigned tmp = _windowDist2_48(begin(f1) + y, begin(f2) + x, parm);
        if (tmp < min)
        {
            min = tmp;
            x_min = x;
        }
    }
    if (min > parm.windowThreshold)
    {
       return 0;
    }
    else 
    {
        if ( x_min - x_pre > parm.med)
        {
            new_cord = _DefaultCord.createCord(create_id_x(genomeId, _DefaultCord.cell2Cord(x_pre + parm.med)),  
                _DefaultCord.cell2Cord(x_pre + parm.med - x_min + y), strand);
        }
        else
        {
            new_cord = _DefaultCord.createCord(create_id_x(genomeId, _DefaultCord.cell2Cord(x_min)), 
                _DefaultCord.cell2Cord(y), strand);
        }
    }
    score = min;
    return new_cord;
}



/*
//WARN::f1 & f2 are supposed to be the same type and parm
uint64_t previousWindow(FeaturesDynamic & f1, 
                        FeaturesDynamic & f2, 
                        uint64_t cord,
                        float & score)
{
    if (f1.isFs1_32())
    {
        return previousWindow1_32(f1.fs1_32, f2.fs1_32, cord, score, *(f1.apx_parm1_32));
    }
    else if (f1.isFs2_48())
    {
        return previousWindow2_48(f1.fs2_48, f2.fs2_48, cord, score, *(f1.apx_parm2_48));
    }
}
uint64_t nextWindow(FeaturesDynamic & f1, 
                    FeaturesDynamic & f2, 
                    uint64_t cord,
                    float & score)
{
    if (f1.isFs1_32())
    {
        return nextWindow1_32(f1.fs1_32, f2.fs1_32, cord, score, *(f1.apx_parm1_32));
    }
    else if (f1.isFs2_48())
    {
        return nextWindow2_48(f1.fs2_48, f2.fs2_48, cord, score, *(f1.apx_parm2_48));
    }
}
*/

uint64_t nextWindow(FeaturesDynamic & f1, 
                    FeaturesDynamic & f2, 
                    uint64_t cord,
                    float & score)
{
    uint64_t genomeId = get_cord_id(cord);
    uint64_t strand = get_cord_strand(cord);
    uint64_t x_pre = _DefaultCord.cord2Cell(get_cord_x(cord));
    uint64_t y_pre = _DefaultCord.cord2Cell(get_cord_y(cord));
    uint64_t x_min = 0;
    uint64_t y;
    uint64_t new_cord = 0;
    unsigned min = ~0;
    ApxMapParmBase * parm;
    unsigned len1, len2; 
    if (f1.isFs2_48())
    {
        parm = f1.apx_parm2_48;
        len1 = length(f1.fs2_48);
        len2 = length(f2.fs2_48);
    }
    else if (f1.isFs1_16())
    {
        parm = f1.apx_parm1_16;
        len1 = length(f1.fs1_16);
        len2 = length(f2.fs1_16);
    }
    else if (f1.isFs1_32())
    {
        parm = f2.apx_parm1_32;
        len1 = length(f1.fs1_32);
        len2 = length(f2.fs1_32);
    }

    if (y_pre + parm->sup * 2 > len1 || x_pre + parm->sup * 2> len2)
        return 0;
    else 
        y = y_pre + parm->med;
    
    for (uint64_t x = x_pre + parm->inf; x < x_pre + parm->sup; x += 1) 
    {
        unsigned tmp = __windowDist(f1, f2, y, x);
        if (tmp < min)
        {
            min = tmp;
            x_min = x;
        }
    }
    if (min > parm->windowThreshold)
    {
       return 0;
    }
    else 
    {
        if ( x_min - x_pre > parm->med)
        {
            new_cord = _DefaultCord.createCord(create_id_x(genomeId, _DefaultCord.cell2Cord(x_pre + parm->med)),  
                _DefaultCord.cell2Cord(x_pre + parm->med - x_min + y), strand);
        }
        else
        {
            new_cord = _DefaultCord.createCord(create_id_x(genomeId, _DefaultCord.cell2Cord(x_min)), 
                _DefaultCord.cell2Cord(y), strand);
        }
    }
    score = min;
    return new_cord;
}

bool extendWindow(FeaturesDynamic & f1, 
                  FeaturesDynamic & f2, 
                  String<uint64_t> & cords, 
                  uint64_t cordy_str,  //extend window between [read_str, read_end) of the read
                  uint64_t cordy_end,  //read_str & read_end are on the same strand of (back(cords))
                  float &  score)
{
    uint cords_p_str = length(cords) - 1;
    uint64_t new_cord = 0;
    while ((new_cord = previousWindow(f1, f2, back(cords), score)) && get_cord_y(new_cord) >= cordy_str)
    {
        appendValue(cords, new_cord);
    } 
    uint cords_p_end = length(cords);
    for (unsigned k = cords_p_str; k < (cords_p_str + cords_p_end) / 2; k++) 
    {
        std::swap(cords[k], cords[length(cords) - k + cords_p_str - 1]);
    }
    while ((new_cord = nextWindow(f1, f2, back(cords), score)) && get_cord_y(new_cord) + window_size < cordy_end)
    {
        appendValue(cords, new_cord);
    }
    return true;    
}

bool initCord(typename Iterator<String<uint64_t> >::Type & it, 
              typename Iterator<String<uint64_t> >::Type hitEnd,
              unsigned & preCordStart,
              String<uint64_t> & cords)
{
    if (empty(cords)){
        appendValue(cords, 0);
        _DefaultHit.setBlockEnd(cords[0]);
    }

    if (it == hitEnd){
        return false;
    }
    else{
        appendValue(cords, *(it));
        //appendValue(cords, _DefaultCord.hit2Cord_dstr(*(it)));
        ++it;
        preCordStart = length(cords) - 1;   
    }
    return true;
}

/*
 * endCord for double strand index
 */
 bool endCord(String<uint64_t> & cord,
              unsigned & preCordStart)
{
    _DefaultCord.setMaxLen(cord, length(cord) - preCordStart);   
    _DefaultHit.setBlockEnd(back(cord));
    return true;
}

/*
 * nextCord for double strand sequence
 * readstr <= nextcody < read_end
 * NOTE::if failed should reteurn new_cord = 0
 */
 uint64_t nextCord(typename Iterator<String<uint64_t> >::Type & it, 
               typename Iterator<String<uint64_t> >::Type const & hitEnd, 
               StringSet<FeaturesDynamic> & f1, 
               StringSet<FeaturesDynamic> & f2,
               unsigned & preCordStart,
               String<uint64_t> & cord,
               uint64_t read_str,   //y bound; on the forward strand
               uint64_t read_end,
               unsigned read_len,
               unsigned thd_cord_size
               )
{
    if (empty(cord))
    {
        appendValue(cord, 0);
        _DefaultHit.setBlockEnd(cord[0]);
    }
    unsigned distThd = getWindowThreshold(f2[0]);
    int f_new_block = 0;
    while (it < hitEnd)
    {
        if (_DefaultHit.isBlockEnd(*(it - 1)))
        {
            _DefaultHit.setBlockEnd(back(cord));
            preCordStart = length(cord);
            f_new_block = 1;
        }
        uint64_t new_cord = *it++;
        if (get_cord_y(new_cord) > get_cord_y(back(cord)) || f_new_block)
        {

            unsigned dist = _windowDist(f1[get_cord_strand(new_cord)], f2[get_cord_id(new_cord)], 
                                        _DefaultCord.cord2Cell(get_cord_y(new_cord)),
                                        _DefaultCord.cord2Cell(get_cord_x(new_cord)));
            uint64_t nyf = (get_cord_strand(new_cord)) ? read_len - 1 - get_cord_y(new_cord) : get_cord_y(new_cord);
            if(dist < distThd && get_cord_y(new_cord) + thd_cord_size < uint64_t(read_len) &&
               nyf >= read_str && nyf + window_size < read_end)
            {
                appendValue(cord, new_cord);
                return new_cord;
            } 
        }
    }
    if (f_new_block)
    {
        _DefaultHit.setBlockEnd(back(cord));
        preCordStart = length(cord);
    }
    return 0;
}

bool path_dst_1(typename Iterator<String<uint64_t> >::Type hitBegin, 
                typename Iterator<String<uint64_t> >::Type hitEnd, 
                StringSet<FeaturesDynamic> & f1,
                StringSet<FeaturesDynamic> & f2, 
                String<uint64_t> & cords,
                uint64_t read_str, //required to be on the forward strand
                uint64_t read_end, // ...forward strand
                uint64_t read_len)
{
    if (hitBegin > hitEnd)
    {
        return false;
    }
    unsigned thd_cord_size = getFeatureWindowSize(f1);
    typename Iterator<String<uint64_t> >::Type it = hitBegin;
    unsigned preBlockPtr;
    float score = 0;
    if(initCord(it, hitEnd, preBlockPtr, cords))
    //if (initCord2(it, hitEnd, f1, f2, preBlockPtr, cords, read_str, read_end, read_len, thd_cord_size))
    {
        do{
            uint64_t strand = get_cord_strand(back(cords));
            uint64_t genomeId = get_cord_id(back(cords));
            uint64_t cordy_str = strand ? read_len - read_end     : read_str;  
            uint64_t cordy_end = strand ? read_len - read_str - 1 : read_end; 
            uint64_t pre_cord_y = _DefaultHit.isBlockEnd(cords[length(cords) - 2]) ?  0 : get_cord_y(cords[length(cords) - 2]) + 1; //left-closed [ , )
            cordy_str = std::max (pre_cord_y, cordy_str);
            extendWindow(f1[strand], f2[genomeId], cords, cordy_str, cordy_end, score);
        }
        while (nextCord(it, hitEnd, f1, f2, preBlockPtr, cords, read_str, read_end, read_len, thd_cord_size));
        set_cord_end (back(cords));
        return endCord(cords, preBlockPtr);   
    }
    set_cord_end (back(cords));
    return false;
}
/*
 * NOTE::cord_str == read_str, while cord_end = read_end - window_size in case of forward strand
 * Do not restrict cordx boundary, since duplications can create cords of cordx[i - 1] > cordx[i]
 */
int path_dst_2(typename Iterator<String<uint64_t> >::Type hitBegin, 
               typename Iterator<String<uint64_t> >::Type hitEnd, 
               StringSet<FeaturesDynamic> & f1,
               StringSet<FeaturesDynamic> & f2, 
               String<uint64_t> & cords,
               uint64_t read_str, //required to be on the forward strand
               uint64_t read_end, // ...forward strand
               uint64_t read_len)
{
    typedef Iterator<String<uint64_t> >::Type Iter; 
    if (hitBegin >= hitEnd - 1) //at leat 2 patterns
    {
        return 0;
    }
    int thd_cord_size = getFeatureWindowSize(f1);
    float score = 0;

    if (empty(cords))
    {
        initCords(cords);
    }

    uint64_t ready_str, ready_end;
    uint64_t cordy_str, cordy_end;
    bool f_sp_l = false; 
    bool f_sp_r = false; 
    bool f_block_str = false;
    bool f_block_end = false;
    bool f_append = false;
    Iterator <String<uint64_t> >::Type itt_next = hitBegin + 1; //hitBegin + 1 < hitEnd for sure 
    Iter itt_first = hitBegin; //point 2 first of block
    for (Iterator<String<uint64_t> >::Type itt = hitBegin; itt < hitEnd; itt = itt_next++) 
    {
        ready_str = get_cord_strand(*itt) ? read_len - read_end : read_str;  
        ready_end = get_cord_strand(*itt) ? read_len - read_str + 1 : read_end; 
        int64_t da_l = isFirstHit(itt) ? 0 : std::abs(int64_t(get_cord_x(*itt) - get_cord_x(*(itt - 1)) - 
            get_cord_y(*itt) + get_cord_y(*(itt - 1))));
        f_sp_l = (da_l > 80) || get_cord_strand((*itt) ^ (*(itt - 1)));
        while (1)
        {   //iterate itt_next within range of thd_cord_size until block_end or sv
            if (itt_next >= hitEnd || isFirstHit(itt_next))
            {
                f_block_end = 1;
                itt_first = itt_next;
                break;
            }
            int64_t da_r = isFirstHit(itt_next) ? 0 : std::abs(int64_t(get_cord_x(*itt_next) - get_cord_x(*(itt_next - 1)) - 
                get_cord_y(*itt_next) + get_cord_y(*(itt_next - 1))));
            f_sp_r = (da_r > 80) || get_cord_strand((*itt_next) ^ (*(itt_next - 1))); //special cord:: add sv conditions here
            if ((get_cord_y(*itt) + thd_cord_size < get_cord_y(*itt_next) && get_cord_x(*itt) + thd_cord_size < get_cord_x(*itt_next)) ||f_sp_r) //out of range or sv 
            {
                break;
            }
            itt_next++;
        }
        if (!f_sp_r && !f_block_end) //normal case 
        {
            cordy_str = f_sp_l ? *itt : (isFirstHit(itt) ? ready_str : get_cord_y(back(cords)));
            cordy_end = get_cord_y(*itt_next);
            appendValue(cords, *itt);
            _DefaultHit.unsetBlockEnd(back(cords));
            f_append = true;
        }
        else
        {
            if (!f_sp_l && get_cord_y(*(itt_next - 1)) >= thd_cord_size && get_cord_x(*(itt_next - 1)) >= thd_cord_size)
            {
                uint64_t new_cord = shift_cord(*(itt_next - 1), -thd_cord_size, -thd_cord_size);
                cordy_str =isFirstHit(itt) ? read_str : get_cord_y(new_cord);
                cordy_end = get_cord_y(*(itt_next - 1)); //don't change it, sv 
                appendValue(cords, new_cord);
                _DefaultHit.unsetBlockEnd(back(cords));
                f_append = true;
            }
            else //skip the pattern
            {
                f_append = false;
            }  
        }
        if (isLastHit(itt) || f_block_end) //if last patttern, force reset cordy_end
        {
            f_block_end = true;
            cordy_end = ready_end;
        }
        if (f_append)
        {
            extendWindow(f1[get_cord_strand(*itt)], f2[get_cord_id(*itt)], cords, cordy_str, cordy_end, score);
            if (f_block_end)
            {
                _DefaultHit.setBlockEnd(back(cords));
            }
        }
        itt_next = f_block_end ? itt_first : itt_next;
        f_sp_l = false;
        f_sp_r = false;
        f_block_end = false;
        f_append = false;
    }
    return 0;
}

/*
 * nextCord for double strand sequence
 * readstr <= nextcody < read_end
 * NOTE::if failed should reteurn new_cord = 0
 */
 uint64_t _filterHits(String<uint64_t> & hits,
                      StringSet<FeaturesDynamic> & f1, 
                      StringSet<FeaturesDynamic> & f2)
{
    unsigned distThd = getWindowThresholdReject(f2);
    int ii_move = 0;
    for (Iterator<String<uint64_t> >::Type it = beginHits(hits); it < endHits(hits); it++)
    {
        unsigned dist = _windowDist(f1[get_cord_strand(*it)], f2[get_cord_id(*it)], 
                                    _DefaultCord.cord2Cell(get_cord_y(*it)),
                                    _DefaultCord.cord2Cell(get_cord_x(*it)));
        if(dist < distThd)
        {
            *(it - ii_move) = *it;
        } 
        else
        {
            ii_move++;
        }
        if (_DefaultHit.isBlockEnd(*it))
        {
            _DefaultHit.setBlockEnd(*(it - ii_move));
        }
    }
    resize (hits, length(hits) - ii_move);
    return 0;
}

bool path_dst(String<uint64_t> & hits,
              StringSet<FeaturesDynamic> & f1,
              StringSet<FeaturesDynamic> & f2, 
              String<uint64_t> & cords,
              uint64_t read_str,
              uint64_t read_end,
              uint64_t read_len,
              int alg_type)
{
    if (isHitsEmpty(hits)) 
    {
        return true;
    }
    if (alg_type == 1){
        return path_dst_1 (beginHits(hits), endHits(hits), f1, f2, cords, read_str, read_end, read_len);
    }
    else if (alg_type == 2){
       _filterHits(hits, f1, f2);
       path_dst_2 (beginHits(hits), endHits(hits), f1, f2, cords, read_str, read_end, read_len);
    }
    return 0;
}

/*==================================================
=          @Section::Mapping and anchoring         =
===================================================*/
/*__________________________________________________
  ---------- @sub::apx chain additionals ----------*/
/*
 * Shortcut of gathering start and end pos for each block of consecutive cords 
 * NOTE::str and end coordniates of each block is shifted by @thd_cord_size / 2 based on the cord coordinates 
 * NOTE::block is different to block of cords.
    1. Block of cords also includes inv cords.
    2. Here block refers to cords of continuous x and y
   inv end will call set_cord_end if @f_set_end is true
 */
int gather_blocks_ (String<uint64_t> & cords, 
                    String<UPair> & str_ends, //result [] closed 
                    String<UPair> & str_ends_p, //result pointer [,) right open
                    uint64_t str_,          //scan [str_, end_) of cords
                    uint64_t end_, 
                    uint64_t read_len,
                    uint64_t thd_large_gap,
                    uint64_t thd_cord_size,
                    int f_set_end, //is set_cord_end for each block 
                    uint64_t (*isEndFunc)(uint64_t),
                    void (*setEndFunc)(uint64_t &)) 
{
    //small gaps processed in the gap module
    clear(str_ends);
    if (length(cords) < 2)
    {
        return 0;
    }
    uint64_t d_shift_max = thd_cord_size / 2;
    uint64_t d_shift = d_shift_max; //NOTE::str end is shifted by this value
    unsigned p_str = str_;
    for (unsigned i = str_ + 1; i < end_; i++)
    {
        if (isEndFunc(cords[i - 1]) || !isCordsConsecutive_(cords[i - 1], cords[i], thd_large_gap))
        {
            d_shift = std::min (read_len - get_cord_y(cords[p_str]) - 1, d_shift_max);
            uint64_t b_str = shift_cord (cords[p_str], d_shift, d_shift);
            d_shift = std::min (read_len - get_cord_y(cords[i - 1]) - 1, d_shift_max);
            uint64_t b_end = shift_cord (cords[i - 1], d_shift, d_shift);
            appendValue (str_ends, UPair(b_str, b_end));
            appendValue (str_ends_p, UPair(p_str, i));
            if (f_set_end)
            {
                setEndFunc(cords[i - 1]);
            }
            p_str = i;
        }
    }
    d_shift = std::min (read_len - get_cord_y(back(cords)) - 1, d_shift_max);
    uint64_t b_str = shift_cord (cords[p_str], d_shift, d_shift);
    d_shift = std::min (read_len - get_cord_y(back(cords)) - 1, d_shift_max);
    uint64_t b_end = shift_cord (back(cords), d_shift, d_shift);
    appendValue (str_ends, UPair(b_str, b_end));
    appendValue (str_ends_p, UPair(p_str, length(cords)));

    return 0; 
}
/*
 * shortcut function
 * Drop too short blocks in cords 
 */
int clean_blocks_ (String<uint64_t> & cords, uint64_t thd_drop_len) 
{
    if (empty(cords)) {return 0;}
    uint64_t ptr = 1;
    uint64_t len = 0;
    for (uint i = 1; i < length(cords); i++)
    {
        len++;
        cords[ptr] = cords[i];
        if (is_cord_block_end (cords[i]))
        {
            ptr = len < thd_drop_len ? ptr - len : ptr;
            len = 0;
        }
        ptr++;
    }
    resize (cords, ptr);
    return 0;
}


/*
 * WARN::1.
   gaps use cord structure to record y coordinate for simplicity. 
   But the x and strand bits are not defined.
   So Do Not use any functions related to x and strand for gaps 
 * 2. Any y of gaps is on the forward strand.
 */
//Collect gaps in coordinates y
int gather_gaps_y_ (String<uint64_t> & cords, 
                    String<UPair> & str_ends,
                    String<UPair> & gaps,
                    uint64_t read_len,
                    uint64_t thd_gap_size)
{
    uint64_t cord_frt = shift_cord(0, 0, 0); //cord at front 
    uint64_t cord_end = shift_cord(0, 0, read_len - 1); //..end

    if (empty(str_ends))
    {
        appendValue (gaps, UPair(cord_frt, cord_end));
        return 0;
    }
    std::sort (begin(str_ends), end(str_ends), [& cords, &read_len](UPair & i, UPair & j)
    {
        uint64_t y1 = get_cord_strand(i.first) ? 
                      read_len - get_cord_y(i.second) - 1 : get_cord_y(i.first);
        uint64_t y2 = get_cord_strand(j.first) ? 
                      read_len - get_cord_y(j.second) - 1 : get_cord_y(j.first);
        return y1 < y2; 
    });

    uint64_t f_cover = 0;
    uint64_t cordy1 = 0;
    uint64_t cordy2 = 0;
    UPair y1 = getUPForwardy (str_ends[0], read_len);
    UPair y2 = y1;
    if (y1.first > thd_gap_size) //check str[0]
    {
        cordy2 = y1.first;
        appendValue(gaps, UPair(cord_frt, cordy2));
    }
    for (unsigned i = 1; i < length(str_ends); i++)
    {
        if (!f_cover)
        {
            y1 = getUPForwardy(str_ends[i - 1], read_len);
            cordy1 = get_cord_y(str_ends[i - 1].second);
        }
        cordy2 = get_cord_y(str_ends[i].first);
        y2 = getUPForwardy(str_ends[i], read_len);
        if (y1.second > y2.second)  
        {
            //skip y2
            //y1.first < y2.first (sorted)
            //=> y2 is all covered by y1;
            f_cover = 1;
        }
        else
        {
            if (y2.first > y1.second &&  //NOTE::uint don't eliminate the first condition
                y2.first - y1.second > thd_gap_size) 
            {
                appendValue (gaps, UPair(cordy1, cordy2));
            }
            f_cover = 0; 
        }
    }
    uint64_t max_y_end = (f_cover) ? y1.second : y2.second;
    if (read_len - max_y_end > thd_gap_size) //be sure y2 = back(str_ends)
    {
        appendValue(gaps, UPair(max_y_end, cord_end));
    }
    return 0;
}

/*
 * Chainable block: y1_end < y2_str, x1_end < x2_str
 * Slightly different on different stands
 * @str_end1.first @str_end1.second required to have same strand
 * @str_end2.first @str_end2.second required to have same strand
 * x are required sorted before call the func: x1_first < x2_first
 */
int _isChainable(UPair str_end1, UPair str_end2,
                 int64_t read_len,
                 int64_t thd_chain_blocks_lower,
                 int64_t thd_chain_blocks_upper)
{
    int64_t dx = get_cord_x(str_end2.first) - get_cord_x(str_end1.second);
    if (get_cord_strand (str_end1.first ^ str_end2.first))
    {
        int64_t dy1 = get_cord_y(str_end2.first) - (read_len - 1 - get_cord_y(str_end1.first));
        int64_t dy2 = read_len - 1 - get_cord_y(str_end2.second) - get_cord_y(str_end1.second);
        return ((dy1 > thd_chain_blocks_lower && dy1 < thd_chain_blocks_upper) ||
                (dy2 > thd_chain_blocks_lower && dy2 < thd_chain_blocks_upper)) &&
                (dx  > thd_chain_blocks_lower && dx < thd_chain_blocks_upper);
    }
    else
    {
        int64_t dy = get_cord_y(str_end2.first - str_end1.second); 
        return dy > thd_chain_blocks_lower && dy < thd_chain_blocks_upper &&
               dx > thd_chain_blocks_lower && dx < thd_chain_blocks_upper;
    }
}
/**
 * sort cords and combine two blocks if they are not overlapped 
 * NOTE::cords of the same block are required to have the same strand
 */
int chainBlocksSimple_ (String<uint64_t> & cords,
                        String<UPair> & str_ends_p,
                        uint64_t read_len,
                        int64_t thd_chain_blocks_lower,
                        int64_t thd_chain_blocks_upper)
{
    if (empty(cords) || empty(str_ends_p))
    {
        return 0;
    }
    std::sort (begin(str_ends_p), end(str_ends_p), [& cords](UPair & i, UPair & j){
        return _DefaultCord.getCordX(cords[i.first]) < _DefaultCord.getCordX(cords[j.first]);
    });
    String<uint64_t> tmp_cords;
    resize (tmp_cords, length(cords));
    unsigned k = 1;
    tmp_cords[0] = cords[0];
    uint64_t cord1 = cords[str_ends_p[0].first];
    uint64_t cord2 = cords[str_ends_p[0].second - 1];
    //UPair y1 = getUPForwardy (UPair(cord1, cord2), read_len);
    //UPair y2 = y1;
    UPair c1 = UPair(cord1, cord2);
    UPair c2 = UPair(cord1, cord2);
    for (unsigned i = 0; i < length(str_ends_p); i++)
    {
        if (i > 0) //skip combine in the first block
        {
            c2 = UPair(cords[str_ends_p[i].first], cords[str_ends_p[i].second - 1]);
            if (_isChainable(c1, c2, read_len, thd_chain_blocks_lower, thd_chain_blocks_upper))
            {
                _DefaultHit.unsetBlockEnd(tmp_cords[k - 1]);
            }
        }
        for (unsigned j = str_ends_p[i].first; j < str_ends_p[i].second; j++)
        {
            tmp_cords[k] = cords[j];
            k++;
        }
        set_cord_end (tmp_cords[k - 1]);
        //y1 = y2;
        c1 = c2;
    }
    cords = tmp_cords;
    return 0;
}

int chainApxCordsBlocks (String<uint64_t> & cords,
                         String<UPair> & str_ends_p,
                         uint64_t read_len,
                         int64_t thd_chain_blocks_lower,
                         int64_t thd_chain_blocks_upper,
                         int alg_type)
{
    if (alg_type == 1)
    {
        chainBlocksSimple_(cords, str_ends_p, read_len, thd_chain_blocks_lower, thd_chain_blocks_upper);
    }
    else if (alg_type == 2)
    {
        //chainBlocksSimple_(cords, str_ends_p, read_len, thd_chain_blocks_lower, thd_chain_blocks_upper);
        int thd_min_chain_len = 1;
        int thd_drop_score = 0;
        ChainScoreMetric chn_score(thd_min_chain_len, thd_drop_score, &getApxChainScore3);
        chainBlocksCords(cords, str_ends_p, chn_score, read_len, 16, 2, &unset_cord_end, &set_cord_end, 1); //2 major chains limit
    }
    return 0;
}

/*________________________________________________
 ----------  streaming index features  ----------*/
unsigned getDIndexMatchAll (DIndex & index,
                            String<Dna5> & read,
                            String<uint64_t> & set,
                            uint64_t read_str,   //map [read_str, read_end) of read
                            uint64_t read_end,
                            MapParms & mapParm)
{
    int dt = 0;
    LShape shape(index.getShape());
    uint64_t xpre = 0;
    hashInit(shape, begin(read));
    for (unsigned k = read_str; k < read_end; k++)
    {
        hashNexth(shape, begin(read) + k);
        uint64_t pre = ~0;
        if (++dt == mapParm.alpha)
        {
            dt = 0;
            if(hashNextX(shape, begin(read) + k) ^ xpre)
            {
                int64_t str_ = queryHsStr(index, shape.XValue);
                int64_t end_ = queryHsEnd(index, shape.XValue);
                //if (end_ - str_ > mapParm.delta || end_ - str_ == 0)
                if (end_ - str_ == 0)
                {
                    continue; 
                }
                for (int64_t i = str_; i < end_; i++)
                {
                    //double t1 = sysTime();
                    uint64_t val = index.getHs()[i];
                    uint64_t cordy = k;
                    uint64_t id = get_cord_id(val);
                    uint64_t strand = FORWARD_STRAND;
                    uint64_t cordx = get_cord_x(val);
                    if (get_cord_strand(val) ^ shape.strand)
                    {
                        cordy = length(read) - 1 - cordy;
                        strand = REVERSE_STRAND;
                    }
                    //NOTE: make sure data structure of anchor is same to cord
                    //!TODO::when val < cordy:  map reads to itself
                    //val - cordy out of bounds.
                    //t1 = sysTime() - t1;
                    //double t2 = sysTime();;
                    if (cordx > cordy)
                    {
                        uint64_t new_anchor = make_anchor(id, cordx, cordy, strand);
                        //todo::may out of boundary
                        //val = shift_cord (val, -cordy, cordy - get_cord_y(val));
                        if (shape2DIndexCordy(shape) == getDIndexCordy(val))
                        {
                            _DefaultHit.setLongPattern(new_anchor);
                        }
                        /*
                        else
                        {
                        }
                        */
                        appendValue(set, new_anchor);
                    }
                    //t2 = sysTime() - t2;
                }
                xpre = shape.XValue;
            }
        }
    }
    return 0;    
}
/**
 * Search double strand pattern in the index and append to anchors
 * NOTE::y value @map_str and @map_end must be coordinates on forward strand of @read rather than its reverse complement
 */
 unsigned getHIndexMatchAll(LIndex & index,
                           String<Dna5> & read,
                           String<uint64_t> & set,
                           uint64_t map_str,
                           uint64_t map_end,
                           MapParms & mapParm)
{   
    int dt = 0;
    LShape shape(index.shape);
    uint64_t xpre = 0;
    hashInit(shape, begin(read));
    uint64_t read_str = get_cord_y(map_str);
    uint64_t read_end = get_cord_y(map_end);
    uint64_t idx_str = _DefaultCord.getCordX(map_str); //NOTE:: id included = id|x
    uint64_t idx_end = _DefaultCord.getCordX(map_end);
    for (unsigned k = read_str; k < read_end; k++)
    {
        hashNexth(shape, begin(read) + k);
        uint64_t pre = ~0;
        if (++dt == mapParm.alpha)
        {
            dt = 0;
            if(hashNextX(shape, begin(read) + k) ^ xpre)
            {
                xpre = shape.XValue;
                uint64_t pos = getXDir(index, shape.XValue, shape.YValue);
                uint64_t ptr = _DefaultHs.getHeadPtr(index.ysa[pos-1]);
                if (index.isEmptyDir(pos) || ptr >= mapParm.delta)
                {
                    dt = 0;
                    continue;
                }
                while ((_DefaultHs.getHsBodyY(index.ysa[pos]) == shape.YValue || 
                        _DefaultHs.getHsBodyY(index.ysa[pos]) == 0))
                {
                    uint64_t idx = _DefaultHs.getHsBodyS(index.ysa[pos]);
                    //if (_DefaultHs.getHsBodyS(pre - index.ysa[pos]) > mapParm.kmerStep &&
                    if(idx >= idx_str && idx < idx_end)
                    {
                        uint64_t id = _getSA_i1(idx);
                        uint64_t x  = _getSA_i2(idx);
                        if (((index.ysa[pos] & _DefaultHsBase.bodyCodeFlag) >>_DefaultHsBase.bodyCodeBit) ^ shape.strand)
                        {
                            uint64_t new_anchor = make_anchor (id, x, length(read) - 1 - k, REVERSE_STRAND);
                            appendValue(set, new_anchor);
                        }
                        else
                        {    
                            uint64_t new_anchor = make_anchor (id, x, k, FORWARD_STRAND);
                            appendValue (set, new_anchor);
                        }
                        pre = index.ysa[pos];
                    }
                    if (++pos > length(index.ysa) - 1)
                    {
                        break;
                    }
                }
            }
        }
    }
    return 0;
}

/*
 * filter anchors in @anchors that can be chained and record i of thoes @anchors[i] in @anchors_list
 */
uint64_t filterAnchorsList(String<uint64_t> & anchors, 
    String<std::pair<unsigned, unsigned> > & anchors_list, 
    uint64_t shape_len, uint64_t thd_anchor_accept_density, uint64_t thd_anchor_accept_min, unsigned thd_anchor_err_bit)
{
    if (length(anchors) <= 1)
    {
        return 0;
    }
    anchors[0] = 0;
    uint64_t thd_1k_bit = 10;
    double t1 = sysTime();
    sort_ska(begin(anchors), end(anchors));
    uint64_t ak2 = anchors[1]; //2/4, 3/4
    uint64_t block_str = 1, sc = 0, count_anchors = 0;
    t1 = sysTime() - t1;
    double t2 = sysTime();
    uint64_t min_y = ULLMAX, max_y = 0;
    for (unsigned i = 1; i < length(anchors); i++)
    {
        uint64_t anc_y = get_cord_y(anchors[i]);
        uint64_t dy2 = std::abs(int64_t(anc_y - get_cord_y(ak2)));
        int f_continuous = (_DefaultCord.getCordX(anchors[i] - ak2) < (dy2 >> thd_anchor_err_bit)); 
        if (f_continuous)
        {
            if (min_y > anc_y) {
                min_y = anc_y;
            }
            if (max_y < anc_y){
                max_y = anc_y;
            }
            ak2 = anchors[(block_str + i) >> 1]; //update the ak to the median 
            ++count_anchors;
        }
        if (!f_continuous || i == length(anchors) - 1)  
        {
            uint64_t thd_accpet_num = std::max(((max_y - min_y) * thd_anchor_accept_density >> thd_1k_bit), thd_anchor_accept_min);
            if (count_anchors > thd_accpet_num)
            {   //anchors[0] is remove
                appendValue(anchors_list, std::pair<unsigned, unsigned>(block_str, i));
            }
            block_str = i;
            ak2 = anchors[i];
            min_y = anc_y;
            max_y = anc_y;
            count_anchors = 1;
        }
    } 
    t2 = sysTime() - t2;
    return 0;
}

/*
 * NOTE::after calling this function,  anchor[0] for additional infos is just removed.
 */
uint64_t filterAnchors1(Anchors & anchors, uint64_t shape_len, uint64_t thd_anchor_accept_density, uint64_t thd_anchor_accept_min, unsigned thd_anchor_err_bit)
{    
    if (anchors.length() <= 1)
    {
        return 0;
    }
    unsigned ii = 0;
    String<std::pair<unsigned, unsigned> > anchors_list;
    double t1 = sysTime();
    filterAnchorsList(anchors.set, anchors_list, shape_len, thd_anchor_accept_density, thd_anchor_accept_min, thd_anchor_err_bit);
    t1 = sysTime() - t1;
    double t2 = sysTime();
    for (unsigned i = 0; i < length(anchors_list); i++)
    {
        for (unsigned j = anchors_list[i].first; j < anchors_list[i].second; j++)
        {
            anchors[ii++] = anchors[j];
        }
    }
    resize (anchors.set, ii);  
    t2 = sysTime() - t2;
    return 0;
}

/*
 * filter by long patterns
 */
uint64_t filterAnchors2(Anchors & anchors, uint64_t shape_len, uint64_t thd_anchor_accept_density, uint64_t thd_anchor_accept_min, unsigned thd_anchor_err_bit, uint64_t thd_max_anchors_num, int64_t thd_anchor_accept_err)
{    
    //double t1 = sysTime();
    if (anchors.length() <= 1)
    {
        return 0;
    }
    String<uint64_t> anchors_long;
    String<std::pair<unsigned, unsigned> > anchors_list;
    String<uint64_t> anchors_xmean;
    resize(anchors_long, anchors.length());
    uint64_t ii = 0;
    for (unsigned i = 0; i < anchors.length(); i++)
    {
        if (_DefaultHit.isLongPattern(anchors[i]))
        {
            anchors_long[ii] = anchors[i];
            _DefaultHit.unsetLongPattern(anchors_long[i]);
            ++ii;
        }
    }
    resize (anchors_long, ii);
    double t1 = sysTime();
    filterAnchorsList(anchors_long, anchors_list, shape_len, thd_anchor_accept_density, thd_anchor_accept_min, thd_anchor_err_bit);
    t1 = sysTime() - t1;
    double t2 = sysTime();
    std::sort (begin(anchors_list), end(anchors_list), [](std::pair<unsigned, unsigned> & a, std::pair<unsigned, unsigned> & b){return a.second - a.first > b.second - b.first;});
    if (length(anchors_list) > thd_max_anchors_num)
    {
        resize (anchors_list, thd_max_anchors_num);
    }
    resize(anchors_xmean, length(anchors_list), 0);

    for (unsigned i = 0; i < length(anchors_list); i++)
    {
        for (unsigned j = anchors_list[i].first; j < anchors_list[i].second; j++)
        {
            anchors_xmean[i] += get_cord_x(anchors_long[j]);
        }
        anchors_xmean[i] /= anchors_list[i].second - anchors_list[i].first;
    }
    t2 = sysTime() - t2;
    double t3 = sysTime();
    ii = 0; //remove anchors[0] 
    for (unsigned i = 0; i < anchors.length(); i++) 
    {
        for (unsigned j = 0; j < length(anchors_xmean); j++) 
        {
            if (std::abs(int64_t(get_cord_x(anchors[i]) - anchors_xmean[j])) < thd_anchor_accept_err)
            {
                anchors[ii] = anchors[i];
                _DefaultHit.unsetLongPattern(anchors[ii]);
                ii++;
                break;
            }
        }
    }
    t3 = sysTime() - t3;
    double t = t1 + t2 + t3;
    resize (anchors.set, ii);
    return 0;
}

uint64_t filterAnchors(Anchors & anchors, uint64_t shape_len, uint64_t thd_anchor_accept_density, uint64_t thd_anchor_accept_min, unsigned thd_anchor_err_bit, uint64_t thd_max_anchors_num, int64_t thd_anchor_accept_err, int alg_type)
{
    if(alg_type == 1)
    {
        filterAnchors1(anchors, shape_len, thd_anchor_accept_density, thd_anchor_accept_min, thd_anchor_err_bit);
    }
    else if (alg_type == 2)
    {
        filterAnchors1(anchors, shape_len, thd_anchor_accept_density, thd_anchor_accept_min, thd_anchor_err_bit);

        //filterAnchors2(anchors, shape_len, thd_anchor_accept_lens, thd_anchor_err_bit, thd_max_anchors_num, thd_anchor_accept_err);
    }
    return 0;
}

uint64_t getDAnchorList(Anchors & anchors, String<int64_t> & list, uint64_t read_str, uint64_t read_end, MapParms & mapParm)
{
    float thd_err_rate = 0.2;
    float thd_anchor_accept_dens = 0.001; //todo::tune err, kmer step related
    float thd_anchor_accept_lens_rate = 0.01;
    int thd_anchor_accept_lens = thd_anchor_accept_lens_rate * (read_end - read_str);
    float thd_anchor_err = 0.2;
    int thd_sig = 10;
    if (anchors.length() <= 1)
    {
        return 0;
    }
    uint64_t ak2 = anchors[0], ak3 = anchors[0]; //2/4, 3/4
    uint64_t c_b = mapParm.shapeLen, sb = 1, sc = 0;
    double t1 = sysTime();
    anchors.sort(anchors.begin(), anchors.end());
    t1 = sysTime() - t1;
    double t2 = sysTime();
    //anchors[0] = anchors[1];
    ak2=anchors[0];
    uint64_t mask = (1ULL << 20) - 1;
    uint64_t min_y = ~0ULL, max_y = 0;

    for (unsigned k = 1; k < anchors.length(); k++)
    {
        int64_t anc_y = get_cord_y(anchors[k]);
        int64_t dy2 = std::abs(int64_t(anc_y - get_cord_y(ak2)));
        int64_t dy3 = std::abs(int64_t(anc_y - get_cord_y(ak3)));
        int f_continuous =  (get_cord_x(anchors[k] - ak2) < thd_anchor_err * dy2 ||
                             get_cord_x(anchors[k] - ak3) < thd_anchor_err * dy3); 
        if (f_continuous)
        {
            int64_t dy = get_cord_y(anchors.set[k]) - get_cord_y(anchors.set[k - 1]);
            dy = std::min(std::abs(dy), int64_t(mapParm.shapeLen));
            c_b += dy; 
            ak2 = anchors[(sb + k) >> 1]; //update the ak to the median 
            ak3 = anchors[k - ((k - sb) >> 2)];
            min_y = std::min(min_y, get_cord_y(anchors[k]));
            max_y = std::max(max_y, get_cord_y(anchors[k]));

        }
        if (!f_continuous || k == anchors.length() - 1)
        {
            if (c_b > thd_anchor_accept_lens && 
                (k - sb) >= uint((max_y - min_y) * thd_anchor_accept_dens))
            {
                anchors.sortPos2(anchors.begin() + sb, anchors.begin() + k);
                appendValue(list, (c_b << 40) + (sb << 20) + k);
            }
            sb = k;
            ak2 = anchors[k];
            ak3 = anchors[k];
            c_b = mapParm.shapeLen;
            min_y = get_cord_y(anchors[k]);
            max_y = get_cord_y(anchors[k]);
        }
    }
    t2 = sysTime() - t2;
    return 0;
}

uint64_t getDHitList(String<uint64_t> & hits, String<int64_t> & list, Anchors & anchors, MapParms & mapParm, int thd_best_n)
{
    uint64_t mask = (1ULL << 20) - 1;
    if (!empty(list))
    {
        std::sort (begin(list), end(list), std::greater<uint64_t>());
        int tmp = length(list) > mapParm.listN ? mapParm.listN : length(list);
        int record_num = 1;
        for (int k = 0; k < tmp; k++)
        {
            if (record_num > thd_best_n)
            {
                break;
            }
            if ((list[0] / 10) < list[k] && list[k])
            {
                unsigned sb = ((list[k] >> 20) & mask);
                unsigned sc = list[k] & mask;
                for (unsigned n = sb; n < sc; n++)
                {
                    appendValue(hits, _DefaultCord.hit2Cord_dstr(anchors[n]));
                    //appendValue(hits, anchors[n]);
                }   
                _DefaultHit.setBlockEnd(back(hits));
                ++record_num;
            }
            else
            {
                break;
            }
        }
        
        return (list[0] >> 40);   
    }
    else
    {
       return 0;
    }
}

uint64_t getDAnchorMatchList(Anchors & anchors, uint64_t read_str, uint64_t read_end, MapParms & mapParm, String<uint64_t> & hit, int thd_best_n)
{
    String <int64_t> list;
    double t1 = sysTime();
    getDAnchorList (anchors, list, read_str, read_end, mapParm);
    t1 = sysTime() - t1;
    double t2 = sysTime();
    getDHitList(hit, list, anchors, mapParm, thd_best_n);
    t2 = sysTime() - t2;
    return 0;
}
/* Remove chains whose y coordinates are coverd by any other longer chain
 * Note::Corresponding str_ends is not required, so update str_ends manually according to the @str_ends_p after calling
 */
int preFilterChains1(String<uint64_t> & hits, String<int>  & hits_score, String<UPair> & str_ends_p)
{
    if (length(str_ends_p) < 2){return 0; }
    //first round remove covers 
    int f_omit = 0;
    int i_move = 0;  //distance to move str_end_p
    int ii_move = 0; //distance to move hits

    for (int i = 1; i < length(str_ends_p); i++)
    {
        f_omit = 0;
        uint64_t stry2 = get_cord_y(hits[str_ends_p[i].first]);
        uint64_t endy2 = get_cord_y(hits[str_ends_p[i].second - 1]);
        for (int j = 0; j < i - i_move; j++)
        {
            uint64_t stry1 = get_cord_y(hits[str_ends_p[j].first]);
            uint64_t endy1 = get_cord_y(hits[str_ends_p[j].second - 1]);
            if (stry1 <= stry2 && endy1 >= endy2) //this block is covered 
            {
                i_move++;
                ii_move += str_ends_p[i].second - str_ends_p[i].first;
                f_omit = 1; 
                break; 
            }
        }
        if (!f_omit && i_move)
        {
            for (int k = str_ends_p[i].first; k < str_ends_p[i].second; k++)
            {
                hits[k - ii_move] = hits[k];
                hits_score[k - ii_move] = hits_score[k];
            }
            str_ends_p[i - i_move].first = str_ends_p[i].first - ii_move;
            str_ends_p[i - i_move].second = str_ends_p[i].second -  ii_move;
        }
    }
    resize (hits, length(hits) - ii_move);
    resize (hits_score, length(hits_score) - ii_move);
    resize (str_ends_p, length(str_ends_p) - i_move);
    return 0;
}

/* Break the chains from the last step into none-overlapperd pieces
 * Note::Corresponding str_ends is not required, so update str_ends manually according to the @str_ends_p after calling
 * getCordXY:: get_cord_y to cut by y coordinates or get_cord_x ...
 */
int preFilterChains2(String<uint64_t> & hits,  String<UPair> & str_ends_p, void (*setEndFunc)(uint64_t &), uint64_t (*getCordXY)(uint64_t))
{
    //second round : break any overlaps
    String<UPair> str_ends_p_tmp;
    String<uint64_t> xy_strs;
    String<uint64_t> xycuts;
    resize (xy_strs, length(str_ends_p));
    resize (xycuts, 2 * length(str_ends_p));
    uint64_t const mask = 1ULL << 62; //a constant large enough
    for (int i = 0; i < length(str_ends_p); i++) //init xycuts
    {
        xycuts[2 * i] = str_ends_p[i].first;
        xycuts[2 * i + 1] = (str_ends_p[i].second - 1) | mask;  //mask to distinguish the .first and .second
    }
    for (int i = 0; i < length(xy_strs); i++) //init
    {
        xy_strs[i] = str_ends_p[i].first;
    }
    std::sort(begin(xycuts), end(xycuts), [&hits, &mask, &getCordXY](uint64_t & a, uint64_t & b){
        return getCordXY(hits[a & (~mask)]) < getCordXY(hits[b & (~mask)]);  //& erase mask sign when sort
    });
    for (int i = 0; i < length(xycuts); i++)
    {
        uint64_t cuty = getCordXY(hits[(xycuts[i] & (~mask))]);
        for (int j = 0; j < length(xy_strs) && xy_strs[j] < length(hits); j++)
        {
            if (cuty < getCordXY(hits[xy_strs[j]])){
                continue;
            }
            for (int k = xy_strs[j]; k < str_ends_p[j].second; k++)
            {
                if (xycuts[i] & mask)
                {
                    if (getCordXY(hits[k]) == cuty)
                    {
                        uint64_t lowery = xy_strs[j];
                        uint64_t uppery = k + 1;
                        if (lowery != uppery)
                        {
                            appendValue(str_ends_p_tmp, UPair(lowery, uppery));
                            xy_strs[j] = uppery;
                        }
                        break;
                    }
                    else if (getCordXY(hits[k]) > cuty)
                    {
                        uint64_t lowery = xy_strs[j];
                        uint64_t uppery = k;
                        if (lowery != uppery)
                        {
                            appendValue (str_ends_p_tmp, UPair(lowery, uppery));
                            xy_strs[j] = uppery;
                        }
                        break;
                    }
                }
                else
                {
                    if (getCordXY(hits[k]) >= cuty)
                    {
                        uint64_t lowery = xy_strs[j];
                        uint64_t uppery = k;
                        if (lowery != uppery)
                        {
                            appendValue (str_ends_p_tmp, UPair(lowery, uppery));
                            xy_strs[j] = uppery;
                        }
                        break;
                    }
                }
            }
        }
    }
    str_ends_p = str_ends_p_tmp;
    std::sort (begin(str_ends_p), end(str_ends_p), [](UPair & a, UPair & b){return a.second < b.second;});
    for (int i = 0; i < length(str_ends_p); i++)
    {
        setEndFunc(hits[str_ends_p[i].second - 1]);
    }

    return 0;
}

//chain anchors to hits
int getAnchorHitsChains(Anchors & anchors, String<uint64_t> & hits, uint64_t shape_len, uint64_t read_len, 
    uint64_t thd_anchor_accept_density, uint64_t thd_anchor_accept_min, uint64_t thd_large_gap, unsigned thd_anchor_err_bit, uint64_t thd_max_anchors_num, int64_t thd_anchor_accept_err, unsigned alg_type_filter) 
{
    double t1 = sysTime();
    //<<debug
    for (int i = 0; i < length(anchors.set); i++)
    {
        uint64_t h = _DefaultCord.hit2Cord_dstr(anchors.set[i]);
        dout << "gahcs" << i << get_cord_y(h) << get_cord_x(h) << "\n";
    }
    //>>debug
    filterAnchors(anchors, shape_len, thd_anchor_accept_density, thd_anchor_accept_min, thd_anchor_err_bit, thd_max_anchors_num, thd_anchor_accept_err, alg_type_filter) ;
    t1 = sysTime() - t1;
    double t2 = sysTime();
    String<UPair> str_ends;
    String<UPair> str_ends_p;
    String<int>   str_ends_p_score;
    String<int> hits_score; 
    initHitsScore(hits_score); //be sure hit_score has the same structure with Hits
    chainAnchorsHits(anchors.set, hits, hits_score);
    gather_blocks_ (hits, str_ends, str_ends_p, 1, length(hits), read_len, thd_large_gap, 0, 0, & is_cord_block_end, & set_cord_end);
    preFilterChains1 (hits, hits_score, str_ends_p);
    preFilterChains2 (hits, str_ends_p, &set_cord_end, &get_cord_y);
    preFilterChains2 (hits, str_ends_p, &set_cord_end, &get_cord_x);
    resize (str_ends_p_score, length(str_ends_p));
    for (int i = 0; i < length(str_ends_p); i++)
    {
        //note the hits_score is in denscending order
        str_ends_p_score[i] = hits_score[str_ends_p[i].first] - hits_score[str_ends_p[i].second - 1];
    }
    chainBlocksHits(hits, str_ends_p, str_ends_p_score, read_len);
    t2 = sysTime() - t2;
    double ts = t1 + t2;
    return 0;
}

/*
 * Map [read_str, read_end) of the read
 */  
uint64_t mnMapReadList(IndexDynamic & index,
                       String<Dna5> & read,
                       Anchors & anchors,
                       uint64_t map_str,
                       uint64_t map_end,
                       MapParms & mapParm,
                       String<uint64_t> & hits,
                       int alg_type,
                       int thd_best_n)
{
    //alg_type = 1;
    uint64_t read_str = get_cord_y(map_str);
    uint64_t read_end = get_cord_y(map_end);
    double tt1 = sysTime();
    if (index.isHIndex())
    {  
        getHIndexMatchAll(index.hindex, read, anchors.set, map_str, map_end, mapParm);    
    }   
    else if (index.isDIndex())
    {
        getDIndexMatchAll(index.dindex, read, anchors.set, read_str, read_end, mapParm);    
    }
    tt1 = sysTime() - tt1;
    double tt2 = sysTime();
    if (alg_type == 1)
    {
        getDAnchorMatchList(anchors, read_str, read_end, mapParm, hits, thd_best_n);
    }
    else if (alg_type == 2)
    {
        double t = sysTime();
        uint64_t thd_anchor_accept_density = 1;// 1 anchor per 1000 bases
        uint64_t thd_anchor_accept_min = 2; //> this
        //uint64_t thd_anchor_accept_lens = (read_end - read_str) * 0.01;
        uint64_t thd_max_anchors_num = 5; //max num of different anchor groups accepted
        int64_t thd_anchor_accept_err = 2500; //(-2500, 2500) 2500=avearge_read_len(10k) * average_err (25%) , use int64_t, not uint64_t
        uint64_t thd_large_gap = 600;     
        unsigned thd_anchor_err_bit = 2; //  = 1/2^2 = 0.25
        int alg_type_filter = 1; 
        if (length(anchors) > 2000)
        {
            alg_type_filter = 2;
        }
        getAnchorHitsChains(anchors, hits, mapParm.shapeLen, length(read), thd_anchor_accept_density, thd_anchor_accept_min, thd_large_gap, thd_anchor_err_bit, thd_max_anchors_num, thd_anchor_accept_err, 2);
        t = sysTime() - t;
        
    }
    tt2 = sysTime() - tt2;
    return 0;
}

/*
 * Approximate map within [@read_str, read_end) of @read
 */
uint64_t apxMap_ (IndexDynamic & index,
                  String<Dna5> & read,
                  Anchors & anchors,
                  MapParms & mapParm,
                  String<uint64_t> & hits, 
                  StringSet<FeaturesDynamic> & f1,
                  StringSet<FeaturesDynamic> & f2,
                  String<uint64_t> & cords, 
                  uint64_t map_str,
                  uint64_t map_end,
                  int alg_type,
                  int thd_best_n) //flag if chain_blocks
{
    clear (hits);
    anchors.init(1);
    initHits(hits);
    mnMapReadList(index, read, anchors, map_str, map_end, mapParm, hits, alg_type, thd_best_n);
    uint64_t read_str = get_cord_y(map_str);
    uint64_t read_end = get_cord_y(map_end);
    path_dst(hits, f1, f2, cords, read_str, read_end, length(read), alg_type);
    return 0;
}

uint64_t apxMap (IndexDynamic & index,
                 String<Dna5> & read,
                 Anchors & anchors,
                 MapParms & mapParm,
                 String<uint64_t> & hit, 
                 StringSet<FeaturesDynamic> & f1,
                 StringSet<FeaturesDynamic> & f2,
                 String<UPair> & apx_gaps,
                 String<uint64_t> & cords, 
                 int f_chain)
{
    double tts = sysTime();
    int64_t thd_cord_size = getFeatureWindowSize(f1); 
    int64_t thd_large_gap = 1000;     // make sure thd_large_gap <= thd_combine_blocks
    int64_t thd_chain_blocks_lower = -100;
    int64_t thd_chain_blocks_upper = 10000; //two blocks of cords will be combined to one if they can be combined and close enough (< this)
    int64_t thd_drop_len = 2;
    thd_drop_len = std::min (thd_drop_len, int64_t(length(read) * 0.05 / thd_cord_size)); //drop blocks length < this
    int thd_best_n = 999; //unlimited best hit;
    clear(apx_gaps);
    if (f_chain)
    {
        MapParms mapParm1 = mapParm;
        MapParms mapParm2 = mapParm;

        mapParm2.alpha = 7;
        mapParm2.listN = mapParm2.listN2;
        thd_best_n = 999; //unlimited best hit;
        int alg_type = 2; //algorithm 1 sort, algorithm 2, dp
        uint64_t map_str = 0ULL;
        uint64_t map_end = create_cord (MAX_CORD_ID, MAX_CORD_X, length(read), 0);
        double t1 = sysTime();
        apxMap_(index, read, anchors, mapParm1, hit, f1, f2, cords, map_str, map_end, alg_type, thd_best_n);
        t1 = sysTime() - t1;
        double t2 = sysTime();
        String<UPair> str_ends;
        String<UPair> str_ends_p;
        clean_blocks_ (cords, thd_drop_len);
        gather_blocks_ (cords, str_ends, str_ends_p, 1, length(cords), length(read), thd_large_gap, thd_cord_size, 1, &is_cord_block_end, &set_cord_end);
        gather_gaps_y_ (cords, str_ends, apx_gaps, length(read), thd_large_gap);
        uint64_t map_d = thd_cord_size >> 1; // cords + to map areas
        uint64_t str_y = 0;                  //stry y of interval between two consecutive blocks
        for (int i = 0; i < length(apx_gaps); i++) //check large gap and re map the gaps
        {
            UPair y = getUPForwardy (apx_gaps[i], length(read));
            uint64_t y1 = y.first;
            uint64_t y2 = y.second;   
            thd_best_n = 1; //best hit only
            map_str =  y1; 
            map_end =  create_cord(MAX_CORD_ID, MAX_CORD_X, y2, 0);
            //apxMap_(index, read, anchors, mapParm2, hit, f1, f2, cords, map_str, map_end, alg_type, thd_best_n);
        }
        
        clear (str_ends);
        clear (str_ends_p);
        gather_blocks_ (cords, str_ends, str_ends_p, 1, length(cords), length(read), thd_large_gap, thd_cord_size, 1, & is_cord_block_end, & set_cord_end);
        t2 = sysTime() - t2;
        double t3 = sysTime();
        chainApxCordsBlocks (cords, str_ends_p, length(read), thd_chain_blocks_lower, thd_chain_blocks_upper, alg_type);
        t3 = sysTime() - t3;
        clean_blocks_ (cords, thd_drop_len);
    }
    else
    {
        float senThr = mapParm.senThr / thd_cord_size;  
        MapParms mapParm1 = mapParm;
        MapParms mapParm2 = mapParm;
        mapParm2.alpha = mapParm2.alpha2;
        mapParm2.listN = mapParm2.listN2;

        int alg_type = 1;
        apxMap_(index, read, anchors, mapParm1, hit, f1, f2, cords, 0, length(read), alg_type, thd_best_n);
        if (_DefaultCord.getMaxLen(cords) < length(read) * senThr)
        {
            clear(cords);
            apxMap_ (index, read, anchors, mapParm2, hit, f1, f2, cords, 0, length(read), alg_type, thd_best_n);
        }   
        clean_blocks_ (cords, thd_drop_len);
    }

    //Mark signs used in the alignment part
    int seg = 0; //add sgn to cords 
    for (int i = 0; i < length(cords); i++)
    {
        set_cord_recd(cords[i], seg);
        set_cord_main(cords[i]);
        if (_DefaultCord.isBlockEnd(cords[i]))
        {
            seg = 1 - seg;
        }
    }
    tts = sysTime() - tts;
    return 0;
}

uint64_t filterGenomes (IndexDynamic & index,
                 String<Dna5> & read,
                 Anchors & anchors,
                 MapParms & mapParm,
                 String<uint64_t> & hit, 
                 StringSet<FeaturesDynamic> & f1,
                 StringSet<FeaturesDynamic> & f2,
                 String<UPair> & apx_gaps,
                 String<uint64_t> & cords, 
                 int f_chain)
{
    int thd_best_n = 999; //unlimited best hit;
    MapParms mapParm1 = mapParm;
    MapParms mapParm2 = mapParm;
    mapParm2.alpha = mapParm2.alpha2;
    mapParm2.listN = mapParm2.listN2;
    int alg_type = 1;
    apxMap_(index, read, anchors, mapParm1, hit, f1, f2, cords, 0, length(read), alg_type,thd_best_n);
    return 0;
}

/*===================================================
=            Features extension for gaps            =
=====================================================*/
bool isOverlap (uint64_t cord1, uint64_t cord2, 
                int revscomp_const, 
                int overlap_size
                )
{
    uint64_t strand1 = _DefaultCord.getCordStrand(cord1);
    uint64_t strand2 = _DefaultCord.getCordStrand(cord2);
    int64_t x1 = _DefaultCord.getCordX(cord1);
    //int64_t y1 = revscomp_const * strand1 - _nStrand(strand1) * get_cord_y(cord1);
    int64_t y1 = _DefaultCord.getCordY (cord1);
    int64_t x2 = _DefaultCord.getCordX(cord2);
    //int64_t y2 = revscomp_const * strand2 - _nStrand(strand2) * get_cord_y(cord2);
    int64_t y2 = _DefaultCord.getCordY (cord2);
    (void) revscomp_const;
    return std::abs(x1 - x2) < overlap_size && 
           std::abs(y1 - y2) < overlap_size && (!(strand1 ^ strand2));
}
/**
 * cord1 is predecessor of cord2 and they are not overlapped
 */
 bool isPreGap (uint64_t cord1, uint64_t cord2, int revscomp_const, int gap_size)
{
    int64_t x1 = _DefaultCord.getCordX(cord1);
    int64_t x2 = _DefaultCord.getCordX(cord2);
    (void) revscomp_const;
    return (x1 + gap_size <= x2);
}
/**
 * cord1 is successor of cord2 and they are not overlapped
 */
 bool isSucGap (uint64_t cord1, uint64_t cord2, int revscomp_const,int gap_size)
{
    return isPreGap (cord2, cord1, revscomp_const, gap_size);
}
/**
 * Extend windows between cord1, and cord2 if they are not overlapped,
   and insert the windows to the k th element of the cords. 
 * If cord1 and cord2 are on the same strand 
   then call previousWindow for cord1 and nextWindow for cord2 until x1 + windowSize < x2
 * If cord1 and cord2 are different strands
   then call nextWindow for cord1 and previousWindow for cord2 along each own strand until it can't be extended any more.
 */         
 int extendPatch(StringSet<FeaturesDynamic> & f1, 
                 StringSet<FeaturesDynamic> & f2, 
                 String<uint64_t> & cords,
                 int kk,
                 uint64_t cord1,
                 uint64_t cord2,
                 int revscomp_const,
                 int overlap_size,
                 int gap_size,
                 uint thd_accept_score
                 )
{
    float score = 0;
    if (isOverlap(cord1, cord2, revscomp_const, overlap_size))
    {
        return 0;
    }
    uint64_t pcord = cord1;
    uint64_t scord = cord2;
    /*
    if (isSucGap(cord1, cord2, revscomp_const))
    {
        pcord = cord2;
        scord = cord1;
    }
    */
    uint64_t strand1 = _DefaultCord.getCordStrand(pcord);
    uint64_t strand2 = _DefaultCord.getCordStrand(scord);
    uint64_t genomeId1 = get_cord_id(pcord);
    uint64_t genomeId2 = get_cord_id(scord);
    int len = 0;
    uint64_t cord = pcord;
    String<uint64_t> tmp;
    uint64_t x_bound = get_cord_x(scord);
    uint64_t y_bound = get_cord_y(scord);
    while (isPreGap(cord, scord, revscomp_const, gap_size))
    {
        cord = nextWindow (f1[strand1], f2[genomeId1], cord, score);
        if (cord && get_cord_y(cord) < y_bound && get_cord_x(cord) < x_bound && score < thd_accept_score)
        {
            appendValue (tmp, cord);
        }
        else
        {
            break;
        }
    }
    uint64_t nw = pcord;
    if (!empty (tmp))
    {
        len += length(tmp);
        nw = back(tmp);
        insert(cords, kk, tmp);
        x_bound = get_cord_x (back(tmp));
        y_bound = get_cord_y (back(tmp));
        clear(tmp);
    }
    else
    {
        x_bound = get_cord_x(pcord);
        y_bound = get_cord_y(pcord);

    }
     
    cord = scord;
    while (isSucGap(cord, nw, revscomp_const, gap_size))
    {
        cord = previousWindow(f1[strand2], f2[genomeId2], cord, score);
        if (cord && get_cord_y (cord) >  y_bound && get_cord_x(cord) > x_bound && score < thd_accept_score)
        {
            appendValue (tmp, cord);
        }
        else
        {
            break;
        }
    }
    if (!empty (tmp))
    {
        for (unsigned i = 0; i < length(tmp) / 2; i++)
        {
            std::swap (tmp[i], tmp[length(tmp) - i - 1]);
        }
        insert(cords, kk + len, tmp);
        len += length(tmp);
    }
    return len;
}


//End all mapper module
