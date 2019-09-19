#include <iostream>
#include "base.h"
#include "cords.h"
#include "shape_extend.h"
#include "index_util.h"
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

//Debug util
int printScript(int64_t & val, CharString header)
{
    int sum = 0;
    std::cout << header << " "; 
    for (int i = 0; i < 16; i++)
    {

        int v = (val >> (i << 2) & 15);
        std::cout << (val >> (i << 2) & 15) << " ";
        sum += v;
    }
    std::cout << " sum=" << sum << "\n";
    return sum;
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
                                unsigned v5,
                                unsigned v6,
                                unsigned v7,
                                unsigned v8) :
    band_width(v1),
    cell_size(v2),
    cell_num(v3),
    windowThreshold(v4),
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
    ApxMapParmBase(0.25, 16, 12, 60, 16, 4, 16, 1000),
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
    ApxMapParmBase(0.25, 16, 12, 36, 16, 4, 32, 1000),
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
    ApxMapParmBase(0.25, 16, 6, 36, 16, 4, 48, 1000)
{
    window_size = windowSize;
}
unsigned getFeatureWindowDelta(FeaturesDynamic & fs)
{
    if (fs.isFs1_32())
    {
        return fs.apx_parm1_32->windowDelta;
    }
    else if (fs.isFs2_48())
    {
        return fs.apx_parm2_48->windowDelta;
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
    if (fs.isFs1_32())
    {
        return fs.apx_parm1_32->windowSize;
    }
    else if (fs.isFs2_48())
    {
        return fs.apx_parm2_48->windowSize;
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
    if (fs.isFs1_32())
    {
        return fs.apx_parm1_32->windowThreshold;
    }
    if (fs.isFs2_48())
    {
        return fs.apx_parm2_48->windowThreshold;
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
    //dout << "parm.scpt_bit\n";// << parm.scpt_bit << parm.scpt_step << parm.scpt_len << "\n";
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
    //dout << "cptn" << parm.scpt_num * parm.scpt_int_step << parm.scpt_int_step << "\n";
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
    /*
    printInt96(*it1, "4841");
    printInt96(*(it1 + 1), "4841");
    printInt96(*(it1 + 2), "4841");
    printInt96(*(it1 + 3), "4841");
    printInt96(*it2, "4842");
    printInt96(*(it2 + 1), "4842");
    printInt96(*(it2 + 2), "4842");
    printInt96(*(it2 + 3), "4842");
    for (int i = 0; i < 12; i++)
    {
        printInt96(*(it1 + i), "int1");
        printInt96(*(it2 + i), "int2");
        dout << "\n";
    }
    return _scriptDist63_31(*it1, *it2) + 
           _scriptDist63_31(*(it1 + 3), *(it2 + 3)) +
           _scriptDist63_31(*(it1 + 6), *(it2 + 6)) + 
           _scriptDist63_31(*(it1 + 9), *(it2 + 9));

    */
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
    //std::cout << "a2i96 " << *it << *(it+1) << " " << ordV << " " << i << " " << addVal << " " << (1 << 255)<< "\n";
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
        //printInt96(buffer[i], "cf3");
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
    //std::cout << "cfs2 " << thd_id << " " << thd_begin << " " << thd_end << "\n";

    int addMod3[3] = {1, 2, 0};   // adMod3[i] = ++ i % 3;
    int96 zero96 = {0, 0, 0};
    std::vector<int96> buffer(3, zero96); //buffer of 3 cells in one script
    setInt96(f[next], zero96);
    //std::cout << "cf48 " << length(f) << "\n";
    for (int i = 0; i < 3; i++) //init buffer and f[next]
    {
        int tmp = thd_begin + (i << parm.scpt_bit); 
        for (int j = tmp; j < tmp + parm.scpt_step; j++)
        {
            add2merInt96(buffer[i], it_str + j);
        }
        //printInt96(buffer[i], "cfp3");
        //std::cout << "cfp3 " << tmp << '\n';
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
            //std::cout << *(it_str + j);
        }
        //printInt96(buffer[ii], "cf2");
        incInt96(f[next], buffer[ii]);
        ii = addMod3[ii];
        next++;
    }
}
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
    score += min;
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
        //<<debug
        //std::cout << "pw1 " << tmp << " " << x * 16 << " " << get_cord_y(cord) << "\n";
        //>>debug
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
    //std::cout << "pw2 " << get_cord_y(new_cord) << " " << min << "\n";
    score += min;
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
    //std::cout << "pw2 " << get_cord_y(new_cord) << " " << min << "\n";
    score += min;
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
    score += min;
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
        //std::cout << "nw1 " << tmp << " " << x * 16 << "\n";
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
    //dout << "nw1" << min << get_cord_y(new_cord) << get_cord_x(new_cord)<< "\n";
    score += min;
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
    score += min;
    std::cout << "nw1 " << get_cord_y(new_cord) << " " << get_cord_x(new_cord) << " " << min << "\n";
    return new_cord;
}

bool extendWindow(FeaturesDynamic & f1, 
                  FeaturesDynamic & f2, 
                  String<uint64_t> & cords, 
                  uint64_t read_str,  //extend window between [read_str, read_end) of the read
                  uint64_t read_end,  //read_str & read_end are on the forward strand 
                  uint64_t read_len,
                  float & score, 
                  uint64_t & strand)
{
    uint64_t pre_cord_y = (_DefaultHit.isBlockEnd(cords[length(cords) - 2]))?
    0:get_cord_y(cords[length(cords) - 2]) + getFeatureWindowDelta(f1);
    unsigned len = length(cords) - 1;
    uint64_t new_cord;
    while (pre_cord_y < get_cord_y(back(cords)))
    {
        new_cord = previousWindow(f1, f2, back(cords), score);
        uint64_t nyf = get_cord_strand(new_cord) ? read_len - 1 - get_cord_y(new_cord) : get_cord_y(new_cord);
        if (new_cord && get_cord_y(new_cord) > pre_cord_y && nyf >= read_str)
        {
            appendValue(cords, new_cord);
        }
        else
        {
            break;
        }
    } 
    for (unsigned k = len; k < ((length(cords) + len) >> 1); k++) 
    {
        std::swap(cords[k], cords[length(cords) - k + len - 1]);
    }
    while (true)
    {
        new_cord = nextWindow(f1, f2, back(cords), score);
        uint64_t nyf = get_cord_strand(new_cord) ? read_len - 1 - get_cord_y(new_cord) : get_cord_y(new_cord);
        if (new_cord && nyf + window_size < read_end)
        {
            appendValue(cords, new_cord);
        }
        else
        {
            break;
        }
    }
    return true;    
}

bool initCord(typename Iterator<String<uint64_t> >::Type & it, 
              typename Iterator<String<uint64_t> >::Type & hitEnd,
              unsigned & preCordStart,
              String<uint64_t> & cord)
{
    if (empty(cord)){
        appendValue(cord, 0);
        _DefaultHit.setBlockEnd(cord[0]);
    }

    if (it == hitEnd){
        return false;
    }
    else{

        appendValue(cord, *(it));
        //appendValue(cord, _DefaultCord.hit2Cord_dstr(*(it)));
        ++it;
        preCordStart = length(cord) - 1;   
    }
    return true;
}

/*
 * endCord for double strand index
 */
 bool endCord(String<uint64_t> & cord,
              unsigned & preCordStart,
              float const & thd_min_block_len)
{
    _DefaultCord.setMaxLen(cord, length(cord) - preCordStart);   
    _DefaultHit.setBlockEnd(back(cord));
    return true;
}

/*
 * nextCord for double strand sequence
 * readstr <= nextcody < read_end
 */
 bool nextCord(typename Iterator<String<uint64_t> >::Type & it, 
               typename Iterator<String<uint64_t> >::Type const & hitEnd, 
               StringSet<FeaturesDynamic> & f1, 
               StringSet<FeaturesDynamic> & f2,
               unsigned & preCordStart,
               String<uint64_t> & cord,
               float const & thd_min_block_len,
               uint64_t read_str,   //y bound; on the forward strand
               uint64_t read_end,
               unsigned readLen,
               unsigned thd_cord_size
               )
{
//TODO: add maxlen of anchor to the first node in cord;
    if (it >= hitEnd)
        return false;
    unsigned distThd;
    unsigned windowDelta = getFeatureWindowDelta(f1);
    while(!_DefaultHit.isBlockEnd(*(it - 1)))
    {
        dout << "gy1" << get_cord_y(*it) << "\n";
        if(get_cord_y(*it) > get_cord_y(back(cord)) +  windowDelta) 
        {
            //uint64_t new_cord = _DefaultCord.hit2Cord_dstr(*it);
            uint64_t new_cord = *it;
            print_cord(new_cord, "gy2");
            uint64_t strand = get_cord_strand(new_cord);
            uint64_t genomeId = get_cord_id(new_cord);
            unsigned dist = _windowDist(f1[strand],                
                                        f2[genomeId], 
                                        _DefaultCord.cord2Cell(get_cord_y(new_cord)),
                                        _DefaultCord.cord2Cell(get_cord_x(new_cord)));
            if (f1[strand].isFs1_32()) 
            {
                distThd = f1[strand].apx_parm1_32->windowThreshold;
            }
            else if (f1[strand].isFs2_48())
            {
                distThd = f1[strand].apx_parm2_48->windowThreshold;
            }
            //dout << "distthd" << dist << distThd << "\n";
            uint64_t nyf = (get_cord_strand(new_cord)) ? readLen - 1 - get_cord_y(new_cord) : get_cord_y(new_cord);
            if(dist < distThd && get_cord_y(new_cord) + thd_cord_size < uint64_t(readLen) &&
               nyf >= read_str && nyf + window_size < read_end)
            {
                appendValue(cord, new_cord);
                ++it;
                return true;
            }
        }
        ++it;
    }
    _DefaultHit.setBlockEnd(back(cord));
    if(it < hitEnd)
    {
        _DefaultCord.setMaxLen(cord, length(cord) - preCordStart);
        preCordStart = length(cord);
        appendValue(cord, *it);
        //appendValue(cord, _DefaultCord.hit2Cord_dstr(*it));
        dout << "gyd222" << back(cord) << "\n";
        print_cord (back(cord), "gyd22");
        ++it;
        return true;
    }
    else
    {
        return false;
    }
}

bool path_dst(typename Iterator<String<uint64_t> >::Type hitBegin, 
              typename Iterator<String<uint64_t> >::Type hitEnd, 
              StringSet<FeaturesDynamic> & f1,
              StringSet<FeaturesDynamic> & f2, 
              String<uint64_t> & cords,
              uint64_t read_str,
              uint64_t read_end,
              uint64_t read_len,
              float const & thd_min_block_len)
{
    unsigned thd_cord_size = getFeatureWindowSize(f1);
    typename Iterator<String<uint64_t> >::Type it = hitBegin;
    unsigned preBlockPtr;
    float score = 0;
    //dout << "win" << _windowDist(begin(f1[0].fs2_48), begin(f2[0].fs2_48)) << "\n";
    if(initCord(it, hitEnd, preBlockPtr, cords))
    {
        do{
            uint64_t strand = get_cord_strand(back(cords));
            uint64_t genomeId = get_cord_id(back(cords));
            extendWindow(f1[strand], f2[genomeId], cords, read_str, read_end, read_len, score, strand);
        }
        while (nextCord(it, hitEnd, f1, f2, preBlockPtr, cords, thd_min_block_len, read_str, read_end, read_len, thd_cord_size));
        set_cord_end (back(cords));
        return endCord(cords, preBlockPtr, thd_min_block_len);   
    }
    set_cord_end (back(cords));
    //std::cout << "[]::path_dist::cord " 
    return false;
}

/*==================================================
=          @Section::Mapping and anchoring         =
===================================================*/

/*__________________________________________________
  ---------- @sub::Generic chain funcs  ----------*/

//Chainning Score metric wrapper: including a score function with corresponding parms.
ChainScoreMetric::ChainScoreMetric(int abort_socre, 
        int(*scoreFunc)(uint64_t const &, uint64_t const &)) 
        : thd_abort_score(abort_socre), getScore(scoreFunc)
        {};
int ChainScoreMetric::getAbortScore()
{
    return thd_abort_score;
}

int getBestChains(String<uint64_t>     & anchors,
                  String<ChainsRecord> & chains,
                  int anchor_end,
                  int (*scoreFunc) (uint64_t const &, uint64_t const &))
{
    if (anchor_end == 0)
    {
        return 0;
    }
    //<<debug
    /*
    dout << "stc<<<<<<<<<<<<<<<<<<<<<<<\n";
    for (auto tmp:anchors)
    {
        uint64_t tmp_cord = _DefaultCord.hit2Cord_dstr(tmp);
        dout << "stc" << get_cord_strand(tmp) <<  get_cord_x(tmp_cord) << "\n";
    }
    */
    //>>debug
    int thd_chain_depth = 20;
    int new_score = 0;
    int new_max_score = 0;
    int max_j = 0;
    int const chain_end_score = 0;
    int const chain_end = -1;
    //ChainsRecord c0 (chain_end_score, 0, chain_end);
    //appendValue(chains, c0);
    chains[0].score = chain_end_score;  
    chains[0].len = 1;
    chains[0].p2anchor = chain_end;
    for (int i = 0; i < anchor_end; i++) 
    {
        int j_str = std::max (0, i - thd_chain_depth);
        max_j = i;
        new_max_score = -1;
        for (int j = j_str; j < i; j++)
        {
            new_score = scoreFunc(anchors[j], anchors[i]);
            //<<debug
                //<< debug
    uint64_t tmp_cord1 = _DefaultCord.hit2Cord_dstr(anchors[j]);
    uint64_t tmp_cord2 = _DefaultCord.hit2Cord_dstr(anchors[i]);
    //dout << "apsss"  << j << i << new_score << get_cord_strand(tmp_cord1) << get_cord_y(tmp_cord1) << get_cord_x(tmp_cord1) << get_cord_strand(tmp_cord2) << get_cord_y(tmp_cord2) << get_cord_x(tmp_cord2) << max_j << new_max_score<< "\n";
    //>>debug
            //>>debug
            if (new_score > 0 && new_score + chains[j].score >= new_max_score)
            {
                max_j = j;
                new_max_score = new_score + chains[j].score;
            }
        }
        if (new_max_score > 0)
        {
            //dout << "apsss_max" << max_j << i << get_cord_y(anchors[i]) << new_max_score << "1507\n";
            chains[i].p2anchor = max_j;
            chains[i].score = new_max_score ;
            chains[i].len = chains[max_j].len + 1;
        }
        else
        {
            chains[i].p2anchor = chain_end;
            chains[i].score = chain_end_score;
            chains[i].len = 1;
        }
    }
    return 0;
}

int createChainsFromAnchors(StringSet<String<uint64_t> > & chains, String<uint64_t> & anchors, int anchor_end, ChainScoreMetric & chn_score1)
{
    if (anchor_end < 2)
    {
        return 0;
    }
    String<uint64_t> chain;
    String<ChainsRecord> chains_record;
    resize (chains_record, anchor_end);
    getBestChains(anchors, chains_record, anchor_end, chn_score1.getScore);
    int delete_score = -1000;
    int bestn = 5;
    int const chain_end = -1;
    for (int i = 0; i < bestn; i++) 
    {
        int max_score = -1;
        int max_str = chain_end;
        int max_len = 0;
        for (int j = 0; j < anchor_end; j++)
        {
            if (chains_record[j].score > max_score)
            {
                max_str   = j;
                max_score= chains_record[j].score;
                max_len   = chains_record[j].len;
            }
        }
        //dout << "maxlen2" << max_len << max_score << max_str << chains_record[max_str].score << "\n";
        if (max_len > 1 && max_score / max_len > chn_score1.getAbortScore())
        {
            clear(chain);
            for (int j = max_str; j != chain_end; j = chains_record[j].p2anchor)
            {
                if (chains_record[j].score != delete_score)
                {
                    //dout << "maxlen2" << get_cord_y(anchors[j]) << "\n";
                    appendValue (chain, anchors[j]);
                    chains_record[j].score = delete_score; 
                }
                else
                {
                    break;
                }
            }
            if (!empty(chain))
            {
                appendValue(chains, chain);
            }
        } 
        if (max_str != chain_end)
        {
            chains_record[max_str].score = delete_score; 
        }
    }   
}

int getApxChainScore(uint64_t const & anchor1, uint64_t const & anchor2)
{
    int64_t dy = get_cord_y(anchor1) - get_cord_y(anchor2);
    if (dy < 0)
    {
        return -10000;
    }
    int64_t thd_min_dy = 50;
    int64_t da = std::abs(int64_t(_DefaultHit.getAnchor(anchor2) - _DefaultHit.getAnchor(anchor1)));
    //<< debug
    uint64_t tmp_cord1 = _DefaultCord.hit2Cord_dstr(anchor1);
    uint64_t tmp_cord2 = _DefaultCord.hit2Cord_dstr(anchor2);
    //dout << "aps" << da << get_cord_strand(tmp_cord1) << get_cord_y(tmp_cord1) << get_cord_x(tmp_cord1) << get_cord_strand(tmp_cord2) << get_cord_y(tmp_cord2) << get_cord_x(tmp_cord2) << "\n";
    //>>debug
    int64_t d_err =  (100 * da) / std::max(dy, thd_min_dy); // 1/100 = 0.01
    //dout << "aps_score x" << dy << d_err << da << get_cord_y(tmp_cord1) << get_cord_y(tmp_cord2) << "\n";
    if (d_err < 10)         {d_err = 0;}
    else if (d_err < 25)    {d_err = 10 + 2 * d_err ;}
    else if (d_err < 100)   {d_err = d_err * d_err / 10 + 40;}
    else                    {d_err = 10000;}

    //d_y
    dy /= 15;
    if (dy < 150)           {dy = dy / 5;}
    else if (dy < 100)      {dy = dy - 30;}
    else if (dy < 10000)    {dy = dy * dy / 200 + 20;}
    else                    {dy = 10000;}
    //dout << "aps_score" << 100 - dy - d_err << dy << d_err << da << get_cord_y(tmp_cord1) << get_cord_y(tmp_cord2) << "\n";
    return 100 - dy - d_err ;    
}

int createApxHitsFromAnchors(String<uint64_t> & hits, String<uint64_t> & anchors)
{
    ChainScoreMetric chn_score(45, &getApxChainScore);
    std::sort(begin(anchors), end(anchors), 
        [](uint64_t & a, uint64_t & b){
            //return _DefaultCord.get_hit_strx(a) > _DefaultCord.get_hit_strx(b);
            return get_cord_x(_DefaultCord.hit2Cord_dstr(a)) > get_cord_x(_DefaultCord.hit2Cord_dstr(b));
        });
    StringSet<String<uint64_t> > chains;
    createChainsFromAnchors(chains, anchors, length(anchors), chn_score); 
    for (auto & chain : chains)
    {
        for (auto & hit : chain) 
        {
            appendValue (hits, _DefaultCord.hit2Cord_dstr(hit));
            //appendValue(hits, hit);
        }
        _DefaultHit.setBlockEnd(back(hits));
    }

    //dout << "anchorlen" << length(anchors) << length(hits) << "\n";
    //print_cords(hits, "hits");
    return 0;
}

void Hit::setBlockStart(uint64_t & val, uint64_t const & flag)
{
    val |= flag;
}

void Hit::setBlockBody(uint64_t & val, uint64_t const & flag)
{
    val &= (~flag);
}

 bool Hit::isBlockStart(uint64_t & val, uint64_t const & flag)
{
    return val & flag;
}

 void Hit::setBlockEnd(uint64_t & val, uint64_t const & flag)
{
    val |= flag;
}

 void Hit::unsetBlockEnd(uint64_t & val, uint64_t const & flag)
{
    val &= ~flag;
}

 void Hit::setBlockStrand(uint64_t & val, uint64_t const & strand, uint64_t const & flag)
{
    if (strand)
        val |= flag;
    else
        val &= ~flag;
}

 bool Hit::isBlockEnd(uint64_t & val, uint64_t const & flag)
{
    return val & flag;
}

 unsigned Hit::getStrand(uint64_t const & val, uint64_t const & flag)
{
    return (val & flag)?1:0;
}

uint64_t Hit::getAnchor(uint64_t const & val)
{
    return (val >> 20 & ((1ULL << 42) - 1) & (~(1ULL << 41)));
}

void _printHit(unsigned j, unsigned id1, unsigned id2, String<uint64_t> & hit, unsigned len)
{
    unsigned end;
    for (unsigned k = 0; k < length(hit); k++)
    {
        if (_DefaultHit.isBlockEnd(hit[k]))
            end = 1;
        else
            end = 0;
        printf("[printhit] %d %d %d %d %d\n", j, id1, id2, len, end);
    }
}

void _printHit(String<uint64_t>  & hit, CharString header)
{
    for (unsigned k = 0; k < length(hit); k++)
    {
        std::cout << "[P]::_printHit() " 
              << get_cord_id(_DefaultCord.hit2Cord(hit[k])) << " " 
              << get_cord_x(_DefaultCord.hit2Cord(hit[k])) << " " 
              << get_cord_y(hit[k]) << "\n";
        if (_DefaultHit.isBlockEnd(hit[k]))
        {
            std::cout << "[P]::_printHit() end\n";

        }
    }
}

unsigned getDIndexMatchAll (DIndex & index,
                            String<Dna5> & read,
                            String<uint64_t> & set,
                            uint64_t read_str,   //map [read_str, read_end) of read
                            uint64_t read_end,
                            MapParm & mapParm)
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
                if (end_ - str_ > mapParm.delta || end_ - str_ == 0)
                {
                    continue; 
                }
                for (int64_t i = str_; i < end_; i++)
                {
                    int64_t val = index.getHs()[i];
                    int64_t cordy = k;
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
                    if (cordx > cordy)
                    {
                        uint64_t new_anchor = make_anchor(id, cordx, cordy, strand);
                        //todo::may out of boundary
                        //val = shift_cord (val, -cordy, cordy - get_cord_y(val));
                        appendValue(set, new_anchor);
                    }
                }
                xpre = shape.XValue;
            }
        }
    }
    //std::cout << "gimad2 " << length(set) << "\n";
    return 0;    
}
/**
 * Search double strand pattern in the index and append to anchors
 * NOTE::@read_str and @read_end must be coordinates of @read rather than its reverse complement
 */
 unsigned getHIndexMatchAll(LIndex & index,
                           String<Dna5> & read,
                           String<uint64_t> & set,
                           uint64_t read_str,
                           uint64_t read_end,
                           MapParm & mapParm)
{   
    int dt = 0;
    LShape shape(index.shape);
    uint64_t xpre = 0;
    hashInit(shape, begin(read));
    for (unsigned k = read_str; k < read_end; k++)
    {
        hashNexth(shape, begin(read) + k);
        uint64_t pre = ~0;
        if (++dt == mapParm.alpha)
//        if (++dt == 1)
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
                    if (_DefaultHs.getHsBodyS(pre - index.ysa[pos]) > mapParm.kmerStep)
                    {
                        uint64_t id = _getSA_i1(_DefaultHs.getHsBodyS(index.ysa[pos]));
                        uint64_t x  = _getSA_i2(_DefaultHs.getHsBodyS(index.ysa[pos]));
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

uint64_t getDAnchorList(Anchors & anchors, String<int64_t> & list, uint64_t read_str, uint64_t read_end, MapParm & mapParm)
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
    anchors.sort(anchors.begin(), anchors.end());
    //anchors[0] = anchors[1];
    ak2=anchors[0];
    uint64_t mask = (1ULL << 20) - 1;
    uint64_t min_y = ~0ULL, max_y = 0;
    //<<debug 
    std::stringstream str;
    str << "anchorlen2 " << anchors.end() - anchors.begin();
    std::cout << str.str()<< "\n";
    //>>debug
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
                //dout << "cbsb" << sb << k << c_b << thd_anchor_accept_lens << uint((max_y - min_y) * thd_anchor_accept_dens) << "\n";
            }
            sb = k;
            ak2 = anchors[k];
            ak3 = anchors[k];
            c_b = mapParm.shapeLen;
            min_y = get_cord_y(anchors[k]);
            max_y = get_cord_y(anchors[k]);
        }

    }
    //print_cords(anchors.set, "set");
}

uint64_t getDHitList(String<uint64_t> & hits, String<int64_t> & list, Anchors & anchors, MapParm & mapParm, int thd_best_n)
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
                //dout << "hits<<<<<<<<<<<\n";
                for (unsigned n = sb; n < sc; n++)
                {
                    appendValue(hits, _DefaultCord.hit2Cord_dstr(anchors[n]));
                    //appendValue(hits, anchors[n]);
                }   
                if (!empty(hits))
                {
                    _DefaultHit.setBlockEnd(back(hits));
                }
                ++record_num;
            }
            else
            {
                break;
            }
        }
        
        //print_cords (hit, "gdialx");
        return (list[0] >> 40);   
    }
    else
    {
       return 0;
    }
}

uint64_t getDAnchorMatchList(Anchors & anchors, uint64_t read_str, uint64_t read_end, MapParm & mapParm, String<uint64_t> & hit, int thd_best_n)
{
    String <int64_t> list;
    getDAnchorList (anchors, list, read_str, read_end, mapParm);
    return getDHitList(hit, list, anchors, mapParm, thd_best_n);
}

/*
 * Map [read_str, read_end) of the read
 */  
uint64_t mnMapReadList(IndexDynamic & index,
                       String<Dna5> & read,
                       Anchors & anchors,
                       uint64_t read_str,
                       uint64_t read_end,
                       MapParm & mapParm,
                       String<uint64_t> & hits,
                       int thd_best_n)
{
    int f_sort = 2;
    if (index.isHIndex())
    {  
        getHIndexMatchAll(index.hindex, read, anchors.set, read_str, read_end, mapParm);    
    }   
    else if (index.isDIndex())
    {
        getDIndexMatchAll(index.dindex, read, anchors.set, read_str, read_end, mapParm);    
    }
    if (f_sort == 1)
    {
        getDAnchorMatchList(anchors, read_str, read_end, mapParm, hits, thd_best_n);
    }
    else if (f_sort == 2)
    {
        erase (anchors.set, 0);
        //<<debug
        for (auto & anchor : anchors.set)
        {
            uint64_t hit = _DefaultCord.hit2Cord_dstr(anchor);
            //dout << "anchors" << get_cord_y(hit) << get_cord_x(hit) << "\n";
        }
        //>>debug
        createApxHitsFromAnchors(hits, anchors.set);
        print_cords(hits, "hits1");
    }
    return 0;
}

/*----------  Chain & Wrapper   ----------*/
/*
 * Shortcut of gathering start and end pos of each block of consecutive cords 
 * NOTE::str and end coordniates of each block is shifted by @thd_cord_size / 2 based on the cord coordinates 
 * NOTE::block is different to block of cords.
    1. Block of cords also includes inv cords.
    2. Here block refers to cords of continuous x and y
   inv end will call set_cord_end if @f_set_end is true
 */
int gather_blocks_ (String<uint64_t> & cords, 
                    String<UPair> & str_ends, //result [] closed 
                    String<UPair> & str_ends_p, //result pointer [,) right open
                    uint64_t readLen,
                    uint64_t thd_large_gap,
                    uint64_t thd_cord_size,
                    int f_set_end) //is set_cord_end for each block 
{
    //small gaps processed in the gap module
    clear(str_ends);
    if (length(cords) < 2)
    {
        return 0;
    }
    uint64_t d_shift_max = thd_cord_size / 2;
    uint64_t d_shift = d_shift_max; //NOTE::str end is shifted by this value
    unsigned p_str = 1;

    for (unsigned i = 2; i < length(cords); i++)
    {
        if (is_cord_block_end(cords[i - 1])||
           !isCordsConsecutive_(cords[i - 1], cords[i], thd_large_gap))
        {
            d_shift = std::min (readLen - get_cord_y(cords[p_str]) - 1, d_shift_max);
            uint64_t b_str = shift_cord (cords[p_str], d_shift, d_shift);
            d_shift = std::min (readLen - get_cord_y(cords[i - 1]) - 1, d_shift_max);
            uint64_t b_end = shift_cord (cords[i - 1], d_shift, d_shift);
            appendValue (str_ends, UPair(b_str, b_end));
            appendValue (str_ends_p, UPair(p_str, i));
            if (f_set_end)
            {
                set_cord_end (cords[i - 1]);
            }
            p_str = i;
        }
    }
    //if (get_cord_y(cords[p_str] - back(cords))) //cords[]_y != back()_y
    //{
        d_shift = std::min (readLen - get_cord_y(back(cords)) - 1, d_shift_max);
        uint64_t b_str = shift_cord (cords[p_str], d_shift, d_shift);
        d_shift = std::min (readLen - get_cord_y(back(cords)) - 1, d_shift_max);
        uint64_t b_end = shift_cord (back(cords), d_shift, d_shift);
        appendValue (str_ends, UPair(b_str, b_end));
        appendValue (str_ends_p, UPair(p_str, length(cords)));
    //}

    return 0; 
}
/*
 * shortcut function
 * Drop too short blocks in cords 
 */
int clean_blocks_ (String<uint64_t> & cords, int64_t thd_drop_len) 
{
    uint p = 1;
    int64_t len = 0;
    for (uint i = 1; i < length(cords); i++)
    {
        len++;
        if (p != i)
        {
            cords[p] = cords[i];
        }
        if (is_cord_block_end (cords[i]))
        {
            if (len < thd_drop_len)
            {
                p -= len;
            }
            len = 0;
        }
        p++;
    }
    resize (cords, p);
}

//shortcut to convert cords pair to y pair (stry, endy) on the forward strand (y projection)
UPair getUPForwardy(UPair str_end, uint64_t readLen)
{
    if (get_cord_strand(str_end.first))
    {
        return UPair(readLen - get_cord_y(str_end.second) - 1,
                     readLen - get_cord_y(str_end.first) - 1);
    }
    else
    {
        return UPair(get_cord_y(str_end.first),
                     get_cord_y(str_end.second));
    }
}
/*
 * WARN::1.
   gaps use cord structure to record y coordinate for simplicity. 
   But the x and strand bits are not defined.
   So Do Not use any functions related to x and strand for gaps 
 * 2. Any y of gaps are are on the forward strand.
 */
//Collect gaps in coordinates y
int gather_gaps_y_ (String<uint64_t> & cords, 
                    String<UPair> & str_ends,
                    String<UPair> & gaps,
                    uint64_t readLen,
                    uint64_t thd_gap_size)
{
    uint64_t cord_frt = shift_cord(0, 0, 0); //cord at front 
    uint64_t cord_end = shift_cord(0, 0, readLen - 1); //..end
    if (empty(str_ends))
    {
        appendValue (gaps, UPair(cord_frt, cord_end));
        return 0;
    }
    std::sort (begin(str_ends), end(str_ends), [& cords, &readLen](UPair & i, UPair & j)
    {
        uint64_t y1 = get_cord_strand(i.first) ? 
                      readLen - get_cord_y(i.second) - 1 : get_cord_y(i.first);
        uint64_t y2 = get_cord_strand(j.first) ? 
                      readLen - get_cord_y(j.second) - 1 : get_cord_y(j.first);
        return y1 < y2; 
    });

    uint64_t f_cover = 0;
    uint64_t cordy1 = 0;
    uint64_t cordy2 = 0;
    UPair y1 = getUPForwardy (str_ends[0], readLen);
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
            y1 = getUPForwardy(str_ends[i - 1], readLen);
            cordy1 = get_cord_y(str_ends[i - 1].second);
        }
        cordy2 = get_cord_y(str_ends[i].first);
        y2 = getUPForwardy(str_ends[i], readLen);
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
    if (readLen - max_y_end > thd_gap_size) //be sure y2 = back(str_ends)
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
                 int64_t readLen,
                 int64_t thd_chain_blocks_lower,
                 int64_t thd_chain_blocks_upper)
{
    int64_t dx = get_cord_x(str_end2.first) - get_cord_x(str_end1.second);
    if (get_cord_strand (str_end1.first ^ str_end2.first))
    {
        int64_t dy1 = get_cord_y(str_end2.first) - (readLen - 1 - get_cord_y(str_end1.first));
        int64_t dy2 = readLen - 1 - get_cord_y(str_end2.second) - get_cord_y(str_end1.second);
        //dout << "_isChainable" << dy1 << dy2 << "\n";
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
int chain_blocks_ (String<uint64_t> & cords,
                   String<UPair> & str_ends_p,
                   uint64_t readLen,
                   int64_t thd_chain_blocks_lower,
                   int64_t thd_chain_blocks_upper)
{
    if (empty(cords) || empty(str_ends_p))
    {
        return 0;
    }
    std::sort (begin(str_ends_p), end(str_ends_p), [& cords](UPair & i, UPair & j){
        return get_cord_x(cords[i.first]) < get_cord_x(cords[j.first]);
    });
    String<uint64_t> tmp_cords;
    resize (tmp_cords, length(cords));
    unsigned k = 1;
    tmp_cords[0] = cords[0];
    uint64_t cord1 = cords[str_ends_p[0].first];
    uint64_t cord2 = cords[str_ends_p[0].second - 1];
    //UPair y1 = getUPForwardy (UPair(cord1, cord2), readLen);
    //UPair y2 = y1;
    UPair c1 = UPair(cord1, cord2);
    UPair c2 = UPair(cord1, cord2);
    for (unsigned i = 0; i < length(str_ends_p); i++)
    {
        if (i > 0) //skip combine in the first block
        {
            //print_cord(cord1, "cb1");
            //print_cord(cord2, "cb2");
            c2 = UPair(cords[str_ends_p[i].first], cords[str_ends_p[i].second - 1]);
            if (_isChainable(c1, c2, readLen, thd_chain_blocks_lower, thd_chain_blocks_upper))
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
    //return 0;
}

/*
 * Approximate map within [@read_str, read_end) of @read
 */
uint64_t apxMap_ (IndexDynamic & index,
                  String<Dna5> & read,
                  Anchors & anchors,
                  MapParm & mapParm,
                  String<uint64_t> & hit, 
                  StringSet<FeaturesDynamic> & f1,
                  StringSet<FeaturesDynamic> & f2,
                  String<uint64_t> & cords, 
                  uint64_t read_str,
                  uint64_t read_end,
                  float cordLenThr,
                  int thd_best_n) //flag if chain_blocks
{
    //todo::wrapper the thds
    clear (hit);
    //clear (anchors);
    anchors.init(1);
    mnMapReadList(index, read, anchors, read_str, read_end, mapParm, hit, thd_best_n);
    path_dst(begin(hit), end(hit), f1, f2, cords, read_str, read_end, length(read), cordLenThr);
}

uint64_t apxMap (IndexDynamic & index,
                 String<Dna5> & read,
                 Anchors & anchors,
                 MapParm & mapParm,
                 String<uint64_t> & hit, 
                 StringSet<FeaturesDynamic> & f1,
                 StringSet<FeaturesDynamic> & f2,
                 String<UPair> & apx_gaps,
                 String<uint64_t> & cords, 
                 float cordLenThr,
                 int f_chain)
{
    int64_t thd_cord_size = getFeatureWindowSize(f1); 
    int64_t thd_large_gap = 1000;     // make sure thd_large_gap <= thd_combine_blocks
    int64_t thd_chain_blocks_lower = -100;
    int64_t thd_chain_blocks_upper = 10000; //two blocks of cords will be combined to one if 1.they can be combined 2. they are close enough (< this)
    int64_t thd_drop_len = 2;
    thd_drop_len = std::min (thd_drop_len, int64_t(length(read) * 0.05 / thd_cord_size)); //drop blocks length < this
    int thd_best_n = 999; //unlimited best hit;
    clear(apx_gaps);
    if (f_chain)
    {
        MapParm mapParm1 = mapParm;
        MapParm mapParm2 = mapParm;
        mapParm2.alpha = 5;
        mapParm2.listN = mapParm2.listN2;
        thd_best_n = 999; //unlimited best hit;
        apxMap_(index, read, anchors, mapParm1, hit, f1, f2, cords, 0, length(read), cordLenThr, thd_best_n);
        String<UPair> str_ends;
        String<UPair> str_ends_p;
        clean_blocks_ (cords, thd_drop_len);
        gather_blocks_ (cords, str_ends, str_ends_p, length(read), thd_large_gap, thd_cord_size, 1);
        gather_gaps_y_ (cords, str_ends, apx_gaps, length(read), thd_large_gap);
        //chain_blocks_ (cords, str_ends, length(read), thd_chain_blocks);
        
        uint64_t map_d = thd_cord_size >> 1; // cords + to map areas
        uint64_t str_y = 0;                  //stry y of interval between two consecutive blocks
        for (int i = 0; i < length(apx_gaps); i++) //check large gap and re map the gaps
        {
            dout << "apxc\n";
            UPair y = getUPForwardy (apx_gaps[i], length(read));
            uint64_t y1 = y.first;
            uint64_t y2 = y.second;
            thd_best_n = 1; //best hit only
            apxMap_(index, read, anchors, mapParm2, hit, f1, f2, cords, y1, y2, cordLenThr, thd_best_n);
        }
        clear (str_ends);
        clear (str_ends_p);
        gather_blocks_ (cords, str_ends, str_ends_p, length(read), thd_large_gap, thd_cord_size, 1);
        chain_blocks_ (cords, str_ends_p, length(read), thd_chain_blocks_lower, thd_chain_blocks_upper);
        clean_blocks_ (cords, thd_drop_len);

    }
    else
    {
        float senThr = mapParm.senThr / thd_cord_size;  
        MapParm mapParm1 = mapParm;
        MapParm mapParm2 = mapParm;
        mapParm2.alpha = mapParm2.alpha2;
        mapParm2.listN = mapParm2.listN2;

        apxMap_(index, read, anchors, mapParm1, hit, f1, f2, cords, 0, length(read), cordLenThr, thd_best_n);
        if (_DefaultCord.getMaxLen(cords) < length(read) * senThr)
        {
            clear(cords);
            apxMap_ (index, read, anchors, mapParm2, hit, f1, f2, cords, 0, length(read), cordLenThr, thd_best_n);
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
    return 0;
}

uint64_t filterGenomes (IndexDynamic & index,
                 String<Dna5> & read,
                 Anchors & anchors,
                 MapParm & mapParm,
                 String<uint64_t> & hit, 
                 StringSet<FeaturesDynamic> & f1,
                 StringSet<FeaturesDynamic> & f2,
                 String<UPair> & apx_gaps,
                 String<uint64_t> & cords, 
                 float cordLenThr,
                 int f_chain)
{
   int thd_best_n = 999; //unlimited best hit;
    MapParm mapParm1 = mapParm;
    MapParm mapParm2 = mapParm;
    mapParm2.alpha = mapParm2.alpha2;
    mapParm2.listN = mapParm2.listN2;
    apxMap_(index, read, anchors, mapParm1, hit, f1, f2, cords, 0, length(read), cordLenThr, thd_best_n);
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
 bool isPreGap (uint64_t cord1, uint64_t cord2, 
                int revscomp_const, 
                int gap_size
                )
{

    int64_t x1 = _DefaultCord.getCordX(cord1);
    int64_t x2 = _DefaultCord.getCordX(cord2);
    (void) revscomp_const;
    return (x1 + gap_size <= x2);
}
/**
 * cord1 is successor of cord2 and they are not overlapped
 */
 bool isSucGap (uint64_t cord1, uint64_t cord2, 
                int revscomp_const,
                int gap_size
                )
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
                 int gap_size
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
        if (cord && get_cord_y(cord) < y_bound && get_cord_x(cord) < x_bound)
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
        if (cord && get_cord_y (cord) >  y_bound && get_cord_x(cord) > x_bound)
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
