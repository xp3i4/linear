#ifndef SEQAN_HEADER_PACMAPPER_H
#define SEQAN_HEADER_PACMAPPER_H

//#include <seqan/seq_io.h>
//#include <seqan/stream.h>
//#include <seqan/index.h>
//#include <seqan/store.h>
//#include <iostream>
//#include <math.h>
//#include <seqan/basic.h>
//#include <bitset>
//#include <climits>

//#include "index_qgram_openaddressing_mn.h"

#include "shape_pm.h"
#include "index_pm.h"
using namespace seqan;

//===================================================================
// const type and variable
//===================================================================
   
//template <typename TSpec = void>
struct Const_{
    
    typedef unsigned _INT_;
    typedef uint64_t _LLT_;
    typedef float _FLT_;
    typedef CharString _CSR_;
    typedef _CSR_   _FILE_;
    typedef _INT_   _SHAPE_;
    typedef seqan::Dna5 _CA5_;

    static const _SHAPE_ _SHAPELEN = 30;
    static const _SHAPE_ _SHAPEWHT = 22;
    static const _INT_ _BLOCKSIZE =  1000;
    static const _INT_ _DELAT = 32; 
    static const _INT_ _THRESHOLD = 30; 
    static const _FLT_ _ALPHA ;
    static const _INT_ _KMERSTEP = 1000;
        
//    void operator (_LLT_ & a, _LLT_ & b)
};

const typename Const_::_FLT_ Const_::_ALPHA = 0.8;

struct Options{
    typename Const_::_INT_  kmerLen;
    typename Const_::_INT_  MiKmLen;
    typename Const_::_FILE_ rFile;
    typename Const_::_FILE_ gFile;
    typename Const_::_FILE_ oFile;
    bool        Sensitive; 
    Options():
        kmerLen(Const_::_SHAPELEN),
        MiKmLen(Const_::_SHAPEWHT),
        rFile(""),
        gFile(""),
        oFile(""),
        Sensitive(false)
        {}
} options;

template <typename TValue>
struct SimpleRecord_
{
    StringSet<typename Const_::_CSR_> ids;
    StringSet<TValue> seq;
    int loadRecord(typename Const_::_FILE_ const & path);
};

template <typename TDna = typename Const_::_CA5_, 
        typename Const_::_SHAPE_ TSpan = Const_::_SHAPELEN,
        typename Const_::_SHAPE_ TWeight = Const_::_SHAPEWHT>
struct PMStruct{
    typedef StringSet<String<TDna> > TStrSet_;
    typedef SimpleRecord_<String<TDna> > TSeq_;
    typedef Minimizer<TSpan, TWeight> Minimizer_;
    typedef Shape<TDna, Minimizer_> TShape_;
    typedef Index<TStrSet_, IndexQGram<Minimizer_, OpenAddressing> >  TIndex_;
};

//template <typename TSpec = void>
struct PMResRecord_{
    typedef typename Const_::_INT_ SeqId_;
    typedef typename Const_::_INT_ SeqLen_;
    typedef typename Const_::_LLT_ MapPos_;
    typedef typename Const_::_LLT_ MapStrand_;
    typedef typename Const_::_LLT_ MapScore_;

    SeqId_      rdId;
    SeqId_      rfId;
    SeqLen_     rdLen;
    SeqLen_     rfLen;
    MapPos_     posStart;
    MapPos_     posEnd;
    MapScore_   score;
    MapStrand_   strand = 1ULL << 63;
};

//template<typename TSpec = void>
struct PMRes{
    typedef String<typename PMResRecord_::SeqId_> ResID;
    //typedef String<typename PMResRecord_::SeqLen_> ResLen;
    typedef String<typename PMResRecord_::MapPos_> ResPos;
    
    ResID   id;
    //ResLen  len;
    ResPos  pos;

    void operator ()(typename PMResRecord_::SeqId_ & id, 
                typename PMResRecord_::MapStrand_ & strand,
                //typename PMResRecord_::SeqanLen & len,
                typename PMResRecord_::MapPos_ & pos);
};

inline void PMRes::operator() (typename PMResRecord_::SeqId_ & mid,
                typename PMResRecord_::MapStrand_ & mstrand,
                //typename PMResRecord_::SeqLen & mlen,
                typename PMResRecord_::MapPos_ & mpos)
{
    appendValue(id, mid);
    appendValue(pos, mpos | mstrand);
}

template <typename TDna = typename Const_::_CA5_, typename TSpec = void>
struct MapperStruct
{
    //typedef PMFile::TReadFile   MReadFile;
    //typedef PMFile::TGnomeFile  MGnomeFile;
    //typedef PMFile::TResultFile MResultFile;

    typedef typename PMStruct<TDna>::TSeq_     MSeq;
    typedef typename PMStruct<TDna>::TSeqSet_  MSeqSet;
    typedef typename PMStruct<TDna>::TShape_   MShape;
    typedef typename PMStruct<TDna>::TIndex_   MIndex;
    typedef typename PMStruct<TDna>::TRes_     MRes;
    
//    typedef PMRes::
};


//template <typename TSpec = void>
struct MapParm{
    typename Const_::_INT_  blockSize;
    typename Const_::_FLT_  alpha;  
    typename Const_::_INT_  delta;
    typename Const_::_INT_  threshold;
    typename Const_::_INT_  kmerStep;
    
    MapParm():
        blockSize(Const_::_BLOCKSIZE),
        alpha(Const_::_ALPHA),
        delta(Const_::_DELAT),
        threshold(Const_::_THRESHOLD),
        kmerStep(Const_::_KMERSTEP)
        {}
    void setMapParm(Options & options);
} _DefaultMapParm;


template <typename TDna = typename Const_::_CA5_, typename TSpec = void>
struct Mapper
{
    typename MapperStruct<TDna>::MSeqSet reads;
    typename MapperStruct<TDna>::MSeqSet gnome;
    typename MapperStruct<TDna>::MShape shape;
    typename MapperStruct<TDna>::MIndex index;

    MapParm        mapParm;
    typename MapperStruct<TDna>::MRes res;
    Mapper()
    {};
    //Mapper(Options const & options)
    //{
    //    loadRecords(options);
    //    setMapParm(options);
    //};
};


int _output()
{
    
    return 0;
}

/*
template <typename Stream, typename Records>
inline void _writeRecords(Stream & stream, Records & records){
      
}

tempalte <typename WOptions_, typename Res>
inline void writeRes(WOptoins & wopt, Res & res)
{
    
}

*/



template <typename TDna>
inline int SimpleRecord_<TDna>::loadRecord(typename Const_::_FILE_ const & path)
{
    //std::cerr << "loading sequence " << std::endl;
    seqan::readRecords(this.ids, this.seq, path);
    return 0; 
}
/*
void loadSequence(Options & options, id)
{
    try{
        if (!_loadSeq())
            throw  
    }
}
*/

/*
struct _compltStr
{
     
    static String<Dna5> _complt = "tgcan";
    void operator()
    {
        resize(res, length(str));
    for (unsigned k = 0; k < length(str); k++)
     //   res[k]=_complt[str[k] - 'A'];
        res[k] = _complt[(unsigned)ordValue(str[k])];

    }
}
*/



const unsigned shapelength = 30; 
const unsigned shapeweight = 22; 
//const unsigned blocklimit = 32;

typedef Iterator<String<Dna> >::Type TIter;
typedef Iterator<String<Dna5> >::Type TIter5;
typedef Shape<Dna, Minimizer<shapelength> > TShape;
typedef Shape<Dna5, Minimizer<shapelength> > TShape5;
typedef Shape<Dna, UngappedShape<shapelength> > TShape_u;
typedef Shape<Dna, SimpleMShape> TMShape;

typedef Index<StringSet<DnaString>, IndexQGram<Minimizer<shapelength>, OpenAddressing > > TIndex;
typedef Index<StringSet<String<Dna5> >, IndexQGram<Minimizer<shapelength>, OpenAddressing > > TIndex5;
typedef Index<StringSet<DnaString>, IndexQGram<UngappedShape<shapelength>, OpenAddressing > > TIndex_u;

static String<Dna5> _complt = "tgcan";
inline void _compltStr(String<Dna5> & str, String<Dna5> & res)
{
    resize(res, length(str));
    for (unsigned k = 0; k < length(str); k++)
     //   res[k]=_complt[str[k] - 'A'];
        res[k] = _complt[(unsigned)ordValue(str[k])];
}

inline void _compltRvseStr(String<Dna5> & str, String<Dna5> & res)
{
    resize(res, length(str));
    for (unsigned k = 0; k < length(str); k++)
     //   res[k]=_complt[str[k] - 'A'];
        res[k] = _complt[(unsigned)ordValue(str[length(str) - k])];
}


template <typename TIndex, typename TObj>
inline unsigned _mnMapReads(TIndex & index, String<TObj> & read,  uint64_t* const anchor, MapParm &  mapParm)
{
    hashInit(index.shape, begin(read));
    for (uint64_t h = 0; h < length(anchor); h++)
        anchor[h] = 0;
    uint64_t x = 1;
    //uint64_t block =(1000 < length(read))?1000:length(read);
    //uint64_t l = block * (0.8/0.2);

    unsigned _block = (mapParm.blockSize < length(read))?mapParm.blockSize:length(read);
    unsigned _dt = _block * (mapParm.alpha / (1 - mapParm.alpha));

    for (unsigned h=0; h <= length(read) - _block; h += _dt)
    {
        hashInit(index.shape, begin(read) + h);
        for (unsigned k = h; k < h + _block; k++)
        {
            hashNext(index.shape, begin(read) + k);
            uint64_t dn = getDir(index, index.shape);
            uint64_t pre = ~0;
            //if(_getBodyCounth(index.dir[dn+1]) - _getBodyCounth(index.dir[dn]) < 32)
            if(_getBodyCounth(index.dir[dn+1]) - _getBodyCounth(index.dir[dn]) < mapParm.delta)
            {
                for (uint64_t n = _getBodyCounth(index.dir[dn]); n < _getBodyCounth(index.dir[dn + 1]); n++)
                {
                    //if (index.sa[n] - pre > 1000)
                    if (index.sa[n] - pre > mapParm.kmerStep)
                    {
                        anchor[x++] = (((index.sa[n])- k) << 20) + k;
                        pre = index.sa[n];
                    }
                }
            }
        }
    }
    anchor[0] = anchor[1];
    uint64_t max = 0, c_b=0, ak=anchor[0], cbb=0, mask_cb = (1<<20) - 1, sb=0;
    unsigned start=0;

    std::sort(begin(anchor), begin(anchor) + x);
    
    for (uint64_t k = 1; k <= x; k++)
    {
        if (anchor[k]-ak < (1000<<20))
        {
            cbb++;
                    }
        else
        {
            std::sort(begin(anchor)+sb, begin(anchor)+k, 
            [&mask_cb](uint64_t &a, uint64_t &b){return (a & mask_cb) < (b & mask_cb);});
            for (uint64_t m = sb+1; m < k; m++)
            {
                if(((anchor[m]-anchor[m-1]) & mask_cb) > shapelength) 
                    c_b += shapelength;
                else
                {
                    c_b += (anchor[m] - anchor[m-1]) & mask_cb; 
                }
            }
            if (c_b > max)
            {
                max = c_b;
                start = sb;
            }
            sb = k;
            ak = anchor[k];
            cbb = 1;
            c_b = shapelength;
        }

    }
    //return anchor[start] >> 20;
    return start + (max << 20) ;
}

template <typename TObj>
void mnMap(StringSet<String<TObj> > & genome, StringSet<String<TObj> > & reads, MapParm & mapParm, String<Pair<uint64_t> > & rs)
{
    //typedef ModifiedString<String<Dna5>, ModComplementDna5> Comp;
    TShape5 shape;
    TIndex5 index(genome);
    String<TObj> comStr;
    double time=sysTime();
    unsigned mask1 = (1<<20)-1, res = 0, tmp = 0;
    uint64_t mask  = 131072 - 1;
    uint64_t anchor[mask + 1] = {1ULL<<63};
    String<uint64_t> result;
    resize(result, length(reads));
    std::cerr << "Creating index \n";
    createQGramIndexDirOnly(index);
    std::cerr << "Filtering reads\n"; 
 //   _compltStr(reads[0], comStr); 
//    std::cout << reads[0] << std::endl;
//    std::cout << comStr << std::endl;
    for (unsigned j = 0; j< length(reads); j++)
    {
        res = _mnMapReads(index, reads[j], anchor, mapParm);

            //std::cout << (res>>20) << std::endl;
        if (res < (mapParm.threshold << 20))
        {
            _compltRvseStr(reads[j], comStr);
            tmp = _mnMapReads(index, comStr, anchor, mapParm);
            res = (tmp > (mapParm.threshold << 20))?tmp:0;
        }
        
        //appendResult();
        
        //std::cout << j << " " << length(reads[j]) << " " << _getSA_i2(anchor[res & mask1] >> 20) << std::endl;
        //appendValue(rs, makePairanchor[res & mask1] >> 20);
        assignValueI1(rs[j], _getSA_i2(anchor[res & mask1] >> 20));
        assignValueI2(rs[j], _getSA_i2(anchor[res & mask1] >> 20) + length(reads[j]));
        
    }
    std::cerr << sysTime() - time << std::endl;
}

//template <typename TSpec = void>
//void printInfo()
//{
//    
//}

template <typename TObj>
String<Pair<uint64_t, uint64_t> > map(StringSet<String<TObj> > & genome, StringSet<String<TObj> > & reads, MapParm & mapParm = _DefaultMapParm)
{
    String<Pair<uint64_t, uint64_t> > result;
    resize(result, length(reads));
    mnMap(genome, reads, mapParm, result);
    for (uint64_t k = 0; k < length(result); k++)
        std::cout << k << " " << result[k].i1 << " " << result[k].i2 << std::endl;
    return result;
}

#endif
