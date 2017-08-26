static const unsigned vec_size = 6;
static const unsigned fShapeLen = 5;
static const unsigned window=16;
static const unsigned delta = window;
static const uint64_t bit=32;


namespace seqan{
struct Feature{
    uint64_t vec[vec_size];
};

const uint64_t m1  = 0x5555555555555555; //binary: 0101...
const uint64_t m2  = 0x3333333333333333; //binary: 00110011..
const uint64_t m4  = 0x0f0f0f0f0f0f0f0f; //binary:  4 zeros,  4 ones ...
const uint64_t m8  = 0x00ff00ff00ff00ff; //binary:  8 zeros,  8 ones ...
const uint64_t m16 = 0x0000ffff0000ffff; //binary: 16 zeros, 16 ones ...
const uint64_t m32 = 0x00000000ffffffff; //binary: 32 zeros, 32 ones
const uint64_t hff = 0xffffffffffffffff; //binary: all ones
const uint64_t h01 = 0x0101010101010101; //the sum of 256 to the power of 0,1,2,3...

inline int popcount64c(uint64_t x)
{
    x -= (x >> 1) & m1;             //put count of each 2 bits into those 2 bits
    x = (x & m2) + ((x >> 2) & m2); //put count of each 4 bits into those 4 bits 
    x = (x + (x >> 4)) & m4;        //put count of each 8 bits into those 8 bits 
    return (x * 0x0101010101010101) >> 56;  //returns left 8 bits of x + (x<<8) + (x<<16) + (x<<24) + ... 
}

inline unsigned _pScore(uint64_t const & x, uint64_t const & y)
{
   return (popcount64c(x ^ y));
}

inline uint64_t getVector(Iterator<String<uint64_t> >::Type const & it)
{
    return (*it << 50 ) + (*(it + 1) << 40 ) + (*(it + 5) << 30) + (*(it + 6) << 20 ) + (*(it + 10) << 10 ) + *(it + 11);
    //return (*it << 50 ) + (*(it + 2) << 40 ) + (*(it + 4) << 30) + (*(it + 6) << 20 ) + (*(it + 8) << 10 ) + *(it + 10);
    //return (*it << 56 ) + (*(it + 2) << 48 ) + (*(it + 4) << 40) + (*(it + 6) << 32 ) + (*(it + 8) << 24 ) + *(it + 16) + *((it) << 8) + *((it));
}

template<typename TIter>
inline void getFeature(TIter sBegin,TIter sEnd, String<uint64_t> & fetrStr)
{
    unsigned count = 0, h = 0;
    uint64_t min = ~0;
    Shape<Dna5,UngappedShape<fShapeLen> > shape;
    hashInit(shape, sBegin);
    for (unsigned j = 0; j < sEnd - sBegin - fShapeLen; j++)
    {
        uint64_t v = hashNext(shape, sBegin + j);
        if (count == window)
        {
            fetrStr[h] = min;   
            min = v;
            count = 0;
            h++;
        }
        if (v < min)
        {
            min = v; 
        }
        count++;
    }
}
/*
template <TIter>
uint64_t _getInt(TIter & const it)
{
    uint64_t t = 0;
    for (unsigned k = 0; k < 12 ; k++)
        t += *(it + k) <<  
}
*/
template <typename TIter>
void feTest(StringSet<String<Dna5> > & reads, TIter const & seqBegin, TIter const & seqEnd)
{
    String<uint64_t> f1, f2;
    std::cout << (seqEnd - seqBegin) /delta << std::endl;
    double time = sysTime();
    uint64_t sum = 0;
    int deltaj = 0;
    resize(f1, (seqEnd - seqBegin) /  delta);
    resize(f2, (seqEnd - seqBegin) / delta);
    for (unsigned h = 2; h < 3; h++)
    //for (unsigned h = 0; h < length(reads); h++)
    {
    getFeature(begin(reads[h]), end(reads[h]), f1);
    getFeature(seqBegin, seqEnd, f2);
    //std::sort(begin(f1), end(f1));
    //std::sort(begin(f2), end(f2));
    std::cout << length(f2) << " " << h << std::endl;
    //unsigned len = (length(f1) < length(f2))?length(f1) - 10:length(f2)-10;
    for (int j = 50; j < length(f2) - 5; j+=10)
    {
        uint64_t vec = getVector(begin(f1) + j);
        uint64_t min =~0; 
        for (int k = j + deltaj - 5; k < j + deltaj + 5; k++)
        { 
            uint64_t tmp = _pScore(vec, getVector(begin(f2) + k));
            if (tmp < min)
            {
                min = tmp;
//                deltaj = k - j;
            }
            //std::cout << " " << _pScore(vec, getVector(begin(f2) + k)) << " ";
            std::cout << " " << tmp << " ";
        } 
    //    std::cout << "\n" << deltaj << "\n" ;
        std::cout << "\n" << std::endl;
    }
    }
    std::cout << sum << std::endl;
    std::cout << sysTime() - time << std::endl;

}

static const float band_width = 0.25;
static const unsigned cmask = ((uint64_t)1<<32) - 1;
static const unsigned cell_size = 16;
static const unsigned cell_num = 12;
static const unsigned window_size = cell_size * cell_num; //16*12
static const unsigned sup = cell_num;
static const unsigned med =ceil((1 - band_width) * cell_num);
static const unsigned inf = ceil((1 - 2 * band_width) * cell_num);

static const unsigned initx = 5; 
static const unsigned inity = 5;

/*
inline uint64_t windowInit(String<uint64_t> &f1, String<uint64_t> & f2,  uint64_t cord)
{
    unsigned min = ~0;
    unsigned minx = 0;
    unsigned miny = 0;
    for (unsigned j = 0; j < inity; j++)
    {
    
        uint64_t v1 = getVector(begin(f1) + j);
    for (unsigned k = cord >> 32; k < (cord >> 32) + initx; k++)
    {
        uint64_t v2 = getVector(begin(f2) + k);
        unsigned tmp = _pScore(v1, v2);
        ////std::cout << tmp << " ";
        if (tmp < min)
        {
            min = tmp;
            minx = k;
            miny = j;
        }
        //std::cout << "window Init " << k << " " << tmp << std::endl;
    }
    }
    return ((uint64_t)minx<< 32) + miny;
}

inline uint64_t windowNext(String<uint64_t> &f1, String<uint64_t> & f2, uint64_t cord)
{
    uint64_t x_pre = cord >> 32;
    uint64_t y_pre = cord & cmask;
    uint64_t x = x_pre + inf;
    uint64_t y = y_pre + med;
    //std::cout << x_pre << " " << y_pre << std::endl; 
    //std::cout << med << " " << inf << std::endl;
    uint64_t v1 = getVector(begin(f1) + y);
    unsigned min = ~0;
    unsigned x_min = 0;
    for (; x < x_pre + sup; x += 1) 
    {
        uint64_t v2 = getVector(begin(f2) + x);
        unsigned tmp = _pScore(v1, v2);
        if (tmp < min)
        {
            min = tmp;
            x_min = x;
        }
    //std::cout << tmp << " " ;
    }
    //std::cout << " window Next " << min << std::endl;
    if ( x_min - x_pre > med)
        return ((uint64_t)(x_pre + med) << 32) + x_pre + med - x_min + y;
    else
        return ((uint64_t)x_min << 32) + y;
}

template <typename TIter>
void pTest(StringSet<String<Dna5> > & reads, TIter const & seqBegin, TIter const & seqEnd)
{
    unsigned b = 33688421;
    //unsigned b = 33688544;
    String<uint64_t> f1, f2;
    std::cout << (seqEnd - seqBegin) /delta << std::endl;
    double time = sysTime();
    int deltaj = 0;
    resize(f1, (seqEnd - seqBegin) /  delta);
    resize(f2, (seqEnd - seqBegin) / delta);
    for (unsigned h = 0; h < length(reads); h++)
    //for (unsigned h = 0; h < 1; h++)
    {
        getFeature(begin(reads[h]), end(reads[h]), f1);
        getFeature(seqBegin, seqEnd, f2);
        //std::cout << h << " window Init " << (windowInit(f1, f2, 0) >> 32) << std::endl;
        //for (unsigned d=0; d < 100; d++)
        //std::cout << d << " "<< d*16+b <<" "<< d*16 << " "<< f1[d] << " " << f2[d] << std::endl;
        //uint64_t cord = windowInit(f1, f2, 0);
        //std::cout <<" " << (cord >> 32) * cell_size + b << " " << ((cord >> 32) + 12)*cell_size  + b << " " << (cord & cmask) * cell_size<< " " << ((cord & cmask) + 12 )* cell_size<< " " << std::endl;

        //for (unsigned j = 0; j < 50; j++)
        //{
        //    cord = windowNext(f1, f2, cord);
        //    std::cout << " " << (cord >> 32) * cell_size + b << " " << ((cord >> 32) + 12)*cell_size  + b << " " << (cord & cmask) * cell_size<< " " << ((cord & cmask) + 12 )* cell_size<< " " << std::endl;
        //} 
            //std::cout << (cord >> 32) << " " << (cord & cmask)<< std::endl;
    }
    //for (unsigned j = 0; j < 50; j++)
    //    std::cout << f2[j] << " ";
    //std::cout << sum << std::endl;
    std::sort(begin(f1), end(f1));
    std::sort(begin(f2), end(f2));
    std::cout << f1[0]*f2[2];
    std::cout << sysTime() - time << std::endl;
}

void pTest1(String<Dna5> genome)
{
    unsigned c[5]={0};
    unsigned l=64;
    for (unsigned j=33688544; j < 33688844;j++) 
    {
        for (unsigned m = 0; m <5; m++)
        {
            //if(j%16==1)
            //{
                std::cout << c[m] << " ";
                if (m==4)
                {
                    //for(unsigned n = 0; n<l;n++)
                    //    std::cout << genome[j+n-1];
                    std::cout <<  j  << std::endl; 
                }
            //} 
            c[m] = 0;
        }
        for (unsigned k=j; k<j+l; k++)
        {
            c[ordValue(genome[k])]++;
        }
    }
}
*/
//======================================================
static const unsigned scriptStep=16;
static const unsigned scriptWindow=6; //2^6
//static const uint64_t scriptMask = (1-3*scriptWindow) -1;
static const int scriptCount[5] = {1, 1<<scriptWindow, 1 <<(scriptWindow * 2), 0, 0};
static const int scriptMask = (1 << scriptWindow) - 1;
static const int scriptMask2 = scriptMask << scriptWindow;



template <typename TIter>
inline uint64_t getScript(TIter const & it)
{
    uint64_t script=0;
     
    for (unsigned k=0; k<(1<<scriptWindow); k++)
    {
        script+= scriptCount[ordValue(*(it+k))];
    }
    return script;
}
/*
inline uint64_t ScriptNext()

template<typename T>
inline _absolute(T t) 
{
    return
}

inline uint64_t deltaScript(uint64_t & s1, uint64_t & s2)
{
    
    return deltaScript();
}
*/


inline int _scriptDist(int const & s1, int const & s2)
{
    return std::abs((s1 & scriptMask)- (s2 & scriptMask)) + std::abs(((s1 & scriptMask2) - (s2 & scriptMask2)) >> scriptWindow) + 
            std::abs((s1>>scriptWindow*2) - (s2>>scriptWindow*2));
}
/*
inline int _scriptMinus(int const & s1, int const & s2)
{
    return ((s1 - s2) & scriptMask) + (((s1>>scriptWindow) - (s2>>scriptWindow)) & scriptMask) +
            ((s1>>(scriptWindow<<1)) - (s2>>(scriptWindow<<1));
}

inline int _scriptPlus(int cont & s1, int const & s2)
{
    return ((s1 + s2) & scriptMask) + (((s1 >> scriptWindow) + (s2 >> scriptWindow)) & scriptMask) 
        + ((s1>>scriptWindow) + (s2>>scriptWindow));
}
*/

template<typename TIter> 
void createFeatures(TIter const & itBegin, TIter const & itEnd, String<int> & f)
{
    unsigned next = 1;
    unsigned window = 1 << scriptWindow;
    f[0] = 0;
    for (unsigned k = 0; k < window; k++)
    {
        f[0] += scriptCount[ordValue(*(itBegin + k))];
    }
    for (unsigned k = scriptStep; k < itEnd - itBegin - window ; k+=scriptStep) 
    {
        f[next] = f[next - 1];
        //for (unsigned j = k; j< k + window; j++)
        //{
        //    std::cout << *(itBegin+j);
        //}
        for (unsigned j = k - scriptStep; j < k; j++)
            f[next] += scriptCount[ordValue(*(itBegin + j + window))] - scriptCount[ordValue(*(itBegin + j))];
        //std::cout << k << " " << (f[next] & scriptMask)<< " " << (f[next] >> 6 & scriptMask) << " " << (f[next] >> 12 & scriptMask) << " "<< (f[next] >> 18) << std::endl;
        next++;
    }

}

template<typename TIter>
inline unsigned _windowDist(TIter const & it1, TIter const & it2)
{
    return _scriptDist(*it1, *it2)+ _scriptDist(*(it1+4),*(it2+4)) + _scriptDist(*(it1+8), *(it2+8));
}

static const unsigned epsilon = 100; 
inline void windowInit(String<int> & f1, String<int> & f2, uint64_t & cord)
{
    unsigned min = ~0;
    uint64_t minx=100;
    std::cout << (cord >> 32) << std::endl;
    for(unsigned k = (cord >>32) - epsilon; k < (cord>>32) + epsilon; k++ ) 
    {
        unsigned tmp = _windowDist(begin(f1), begin(f2) + k);
        //std::cout << k << " " << tmp << " " << min << std::endl;
        if (tmp < min)
        {
            min = tmp;
            minx = k;
        }
    }
    std::cout << "Init min " << min << std::endl;
    cord = (cord & cmask) + (minx << 32);
}

inline uint64_t windowNext(String<int> &f1, String<int> & f2, uint64_t & cord)
{
    uint64_t x_pre = cord >> 32;
    uint64_t y_pre = cord & cmask;
    uint64_t y = y_pre + med;
    unsigned min = ~0;
    unsigned x_min = 0;
    for (uint64_t x = x_pre + inf; x < x_pre + sup; x += 1) 
    {
        unsigned tmp = _windowDist(begin(f1) + y, begin(f2) + x);
        //std::cout << x << " " << tmp << std::endl;
        if (tmp < min)
        {
            min = tmp;
            x_min = x;
        }
    }
        //std::cout << "min= " << min << std::endl;
    if ( x_min - x_pre > med)
        cord = ((uint64_t)(x_pre + med) << 32) + x_pre + med - x_min + y;
    else
        cord = ((uint64_t)x_min << 32) + y;
    return min;
}

/*
template <tyepname TIter>
void getPath(TIter const & it1Begin, TIter const & it1End, TIter const & it2Begin, TIter const &it2End)
{
     
}
*/

void pTest4(StringSet<String<Dna5> > & reads, StringSet<String<Dna5> > & genomes, StringSet<CharString> & id_r, uint64_t a[])
{
    std::cout << "done \n" << length(a);
    unsigned window = 1 << scriptWindow;
    String<int> f1, f2;
    double time = sysTime();
    resize(f2, ((end(genomes[0]) -begin(genomes[0]) - window) / scriptStep));
    createFeatures(begin(genomes[0]), end(genomes[0]), f2);
    uint64_t cord = 20ULL << 32;
    std::cout << sysTime() - time<< std::endl;
    for (unsigned j=0; j < length(a); j++) 
    {
        resize(f1, (length(reads[j]) - window) / scriptStep);
        createFeatures(begin(reads[j]), end(reads[j]), f1);
        //cord = (uint64_t)a[j] << 32;
        cord = a[j];
//        windowInit(f1, f2, cord); 
        
        //for (unsigned k = 0; k < length(reads[j]); k++)
        //    std::cout << genomes[0][a[j] * 16 + k];
        //std::cout << j << std::endl; 
        //std::cout << a[j] * 16 << " " << (cord >> 32) << " " << (cord & cmask) << std::endl;
        std::cout << id_r[j] << std::endl;
        while ((cord &cmask) < (length(reads[j]) - window) / scriptStep - 10)
        {
            
        std::cout << " min =  " << (windowNext(f1, f2, cord)) ;//<< "   " << (cord >> 32) * 16 << " " << (cord & cmask) * 16<< std::endl ;
        }
        //std::cout << std::endl;
    }
        //std::cout << (cord >> 32) << " " << (cord & cmask) << std::endl;
    std::cout <<  std::endl << sysTime() - time << std::endl;
}

void dTest(StringSet<String<Dna5> > & reads, StringSet<String<Dna5> > & genome, String<uint64_t> const & a) //test script distribution
{
    unsigned window = 1 << scriptWindow;
    String<int> f1, f2;
    double time = sysTime();
    uint64_t amask = (1ULL << 32) - 1;
    resize(f2, ((end(genome[0]) -begin(genome[0]) - window) / scriptStep));
    createFeatures(begin(genome[0]), end(genome[0]), f2);    
    unsigned count = 0;
    for (unsigned j = 0; j < length(a); j++)
    {
        //std::cout << j << " " << (a[j] >> 32)  * 16 << " " << (a[j] & amask) * 16<< " ";
        if (a[j]==0)
        {
            std::cout << " J " << j << std::endl;
            resize(f1, (length(reads[count]) - window) / scriptStep);
            createFeatures(begin(reads[count]), end(reads[count]), f1);
            count += 1;
        }
        else 
            std::cout << (a[j] >> 32) * 16<< " " << (a[j] & amask) * 16 << " " << _windowDist(begin(f1) + (a[j]>>32), begin(f2) + (a[j] & amask)) << std::endl;
    }
}

template <typename TIter>
void pTest3(StringSet<String<Dna5> > & reads, StringSet<String<Dna5> > & genomes, TIter const & itBegin, TIter const & itEnd)
{
    unsigned window = 1 << scriptWindow;
    String<int> f1, f2;
    double time = sysTime();
    resize(f2, ((itEnd - itBegin - window) / scriptStep));
    createFeatures(itBegin, itEnd, f2);
    uint64_t cord = 20ULL << 32;
    for (unsigned j=0; j < length(reads); j++) 
    {
        resize(f1, (length(reads[j]) - window) / scriptStep);
        createFeatures(begin(reads[j]), end(reads[j]), f1);
        windowInit(f1, f2, cord); 
        std::cout << (cord >> 32) << " " << (cord & cmask) << std::endl;
        while ((cord &cmask) < (length(reads[j]) - window) / scriptStep)
        {
            (windowNext(f1, f2, cord));
//        std::cout << (cord >> 32) << " " << (cord & cmask) << std::endl;
        }
    }
        std::cout << (cord >> 32) << " " << (cord & cmask) << std::endl;
    std::cout <<  std::endl << sysTime() - time << std::endl;
}

template <typename TIter>
void pTest2(StringSet<String<Dna5> > & reads, TIter const & itBegin, TIter const & itEnd)
{
    unsigned window = 1 << scriptWindow;
    String<int> f1, f2;
    resize(f2, ((itEnd - itBegin - window) / scriptStep));
    createFeatures(itBegin, itEnd, f2);
    int t1 = scriptCount[0]*20 + scriptCount[1]*13 + scriptCount[2]*17;
    int t2 = scriptCount[0]*19 + scriptCount[1]*15 + scriptCount[2]*13;
    std::cout << std::endl;
    double time = sysTime();
    for (unsigned j=0; j < length(reads); j++) 
    {
        resize(f1, (length(reads[j]) - window) / scriptStep);
        createFeatures(begin(reads[j]), end(reads[j]), f1);
 
        for (unsigned k=0; k<length(reads[j]) - window; k+=150)
        {
            int s = 0;
            for (unsigned h = k-5; h<k+5; h++)
            { 
                s =  _scriptDist(f1[k], f2[h])+ _scriptDist(f1[k+4],f2[h+4]) + _scriptDist(f1[k+8], f2[h+8]);
            }
            
          std::cout << k << " " << s << std::endl;
        }
    }
    std::cout <<  std::endl << sysTime() - time << std::endl;
}


}
