#include "base.h"

using namespace seqan;

//===================================================================
// variable and type def
//===================================================================

const float    base_alpha_      = 0.75;
const unsigned base_shape_len_  = 25;
const unsigned base_block_size_ = 100;
const unsigned base_delta_      = 32; 
const unsigned base_threshold_  = 30; 
const unsigned base_kmer_step_  = 1000;
const uint64_t base_llt_max_    = ~0;

using std::cerr;

Options::Options():
        oPath(""),
        sensitivity(1),
        thread(16),
        index_t(1),
        feature_t(2),
        gap_len(0),
        aln_flag(0)
        {
           date = __DATE__; 
        }
std::string Options::getOutputPath() const {return oPath;};

/*
 * flip strand from 0, 1 to -1, 1;
 * strand = 0, 1, other values is not allowed
 * return -1 , 1
 */
uint64_t _nStrand(uint64_t strand)
{
    return (strand << 1) - 1;
}
/*
 * flip the coordinates if strand = -1
 * do nothing if strand = 0
 * len is the length of the sequences;
 */
uint64_t _flipCoord (uint64_t coord, uint64_t len, uint64_t strand)
{
    return len * strand - _nStrand(strand) * coord;
}

MapParm::MapParm():
        blockSize(base_block_size_),
        delta(base_delta_),
        threshold(base_threshold_),
        kmerStep(base_kmer_step_),
        shapeLen(base_shape_len_),
        sensitivity(0),
        anchorDeltaThr(),
        minReadLen(1000),
        listN(1),
        listN2(1),
        alpha(base_alpha_),
        alpha2(0.5),
        anchorLenThr(0.02),                  // anchors with lenghth > this parameter is pushed into the queue
        rcThr(0.8),                        // when max anchors in the queue with length < this parameters, reverse complement search will be conducted
        cordThr(0.8),
        senThr(0.8),
        clsThr(0.1)
    {}
MapParm::MapParm(unsigned bs, unsigned dt, unsigned thr, 
            unsigned ks, unsigned sl, unsigned st,
            unsigned ad, unsigned mr, unsigned listn,
            unsigned listn2,
            float ap, float ap2, float alt, float rt, float ct, float sent, float clst):
        blockSize(bs),
        delta(dt),
        threshold(thr),
        kmerStep(ks),
        shapeLen(sl),
        sensitivity(st),
        anchorDeltaThr(ad),
        minReadLen(mr),
        listN(listn),
        listN2(listn2),
        alpha(ap),
        alpha2(ap2),
        anchorLenThr(alt),                  // anchors with lenghth > this parameter is pushed into the queue
        rcThr(rt),                        // when max anchors in the queue with length < this parameters, reverse complement search will be conducted
        cordThr(ct),
        senThr(sent),
        clsThr(clst)
        {} 
MapParm::MapParm(MapParm & parm):
        blockSize(parm.blockSize),
        delta(parm.delta),
        threshold(parm.threshold),
        kmerStep(parm.kmerStep),
        shapeLen(parm.shapeLen),
        sensitivity(parm.sensitivity),
        anchorDeltaThr(),
        minReadLen(parm.minReadLen),
        listN(parm.listN),
        listN2(parm.listN2),
        alpha(parm.alpha),
        alpha2(parm.alpha2),
        anchorLenThr(parm.anchorLenThr),
        rcThr(parm.rcThr),
        cordThr(parm.cordThr),
        senThr(parm.senThr),
        clsThr(parm.clsThr)
        {}

std::ifstream::pos_type _filesize(const char* filename)
{
    std::ifstream in(filename,  std::ifstream::binary | std::ifstream::ate);
    return in.tellg(); 
}

int readRecords_block (StringSet<CharString> & ids, StringSet<String<Dna5> > & reads, String<size_t> & lens, SeqFileIn & fin, int blockSize)
{
    int start = length(reads);
    readRecords(ids, reads, fin, blockSize);
    for (unsigned k = 0; k < length(reads) - start; k++)
    {
        appendValue (lens, length(reads[start + k]));
    }
    return 0;
}

/*
 *[]::load all records in one genome file
 */
std::pair<uint, uint> loadRecords(StringSet<String<Dna5> > & seqs, 
            StringSet<CharString> & ids, 
            Options::PathType path
            )
{
    double time = sysTime();
    SeqFileIn gFile;
    std::pair<uint, uint> res;
    if (!open(gFile, toCString(path)))
    {
        serr.print_message("\033[1;31mError:\033[0m can't open file ", 2, 0, std::cerr);
        serr.print_message(toCString(path), 0, 1, std::cerr);
        res =std::make_pair (uint(~0), uint(~0));
        return res;
    }
    double fileSize = _filesize (toCString(path));
    bool flag = false;
    unsigned seqCount = 0;
    double currentFileSize = 0;
    StringSet<String<char> > dotstatus;
    resize(dotstatus, 3);
    dotstatus[0] = ".   ";
    dotstatus[1] = "..  ";
    dotstatus[2] = "... ";
    unsigned len_sum = 0;
    int error = 0;
#pragma omp parallel
{
    #pragma omp sections
    {
        #pragma omp section
        {
            //unsigned preSeqCount = 0;
            String <char> probar;
            float prepercent = 0, percent = 0, showpercent = 0, v = 0.87 ;
            unsigned k = 1;
            while (!flag)
            {
                prepercent = percent;
                percent = currentFileSize / fileSize * 100;
                percent = (percent > 100)?prepercent:percent;
                showpercent += v;
                showpercent = (showpercent > percent)?percent:showpercent;
                
                std::cerr << "                                                            \r";
                if (seqCount > 2)
                {
                    std::cerr << "=>Read genomes" << dotstatus[(k - 1)/10 %3] << "            " << seqCount << "/" << std::setprecision(2) << std::fixed << showpercent << "%\r";
                }
                else
                {
                    std::cerr << "=>Read genomes" << dotstatus[(k - 1)/10 %3] << "\r";
                }
                std::this_thread::sleep_for(std::chrono::milliseconds(100));
                k++;
            }
        }
        #pragma omp section
        {
            PMRecord::RecId tmp_id;
            PMRecord::RecSeq tmp_seq;
            while (!atEnd(gFile))
            {
                clear (tmp_id);
                clear (tmp_seq);
                try
                {
                    readRecord (tmp_id, tmp_seq, gFile);
                }
                catch (Exception const & e)
                {
                    std::string msg1 = "File: " + path + " ";
                    serr.print_message (msg1, 2, 0, std::cerr);
                    serr.print_message ("[", 20, 0, std::cerr);
                    serr.print_message("\033[1;31mError:\033[0m can't read records in file]", 0, 1, std::cerr);
                    error = 1;
                }
                currentFileSize += length(tmp_seq);
                appendValue (ids, tmp_id);
                appendValue (seqs, tmp_seq);
                len_sum += length(tmp_seq);
                ++seqCount;
            }
            flag = true;
        }
    }
}
    if (error)
    {
        res = std::make_pair(uint(~0), uint(~0));
    }
    else
    {
        res = std::make_pair (len_sum, seqCount);
    }
    return res;
}

int loadRecords(StringSet<String<Dna5> > & seqs, 
                StringSet<CharString> & ids, 
                Options::PathsType & paths)
{
    int status = 0;
    for (uint i = 0 ; i < length(paths); i++)
    {
        double time = sysTime();
        std::pair<uint, uint> res = loadRecords(seqs, ids, paths[i]);
        uint len_sum = res.first;
        uint seqCount = res.second; 
        if (len_sum == ~0 || seqCount == ~0)
        {
            status += 1;
            continue;
        }
        if (i == 0)
        {
            serr.print_message ("--Read genomes      ", 0, 1, cerr);
        }
        std::string msg1 = "File: " + paths[i] + " ";
        serr.print_message (msg1, 2, 0, cerr);

        serr.print_message ("[", 20, 0, cerr);
        serr.print_message (double(seqCount), 0, 0, cerr);
        serr.print_message (" sequences; ", 0, 0, cerr);

        serr.print_message (double(len_sum >> 20), 0, 0, cerr);
        serr.print_message (" mbases; ", 0, 0, cerr);

        std::string msg3 = "Elapsed time[s] ";
        serr.print_message (msg3, 0, 0, cerr);
        serr.print_message (sysTime() - time, 0, 0, cerr);
        
        serr.print_message ("\033[1;32m 100%\033[0m", 0, 0, cerr);

        serr.print_message ("]", 0, 1, cerr);
        //serr.print_message ()

    }
    return status;
}

PMRecord::PMRecord(Options & options)
{

}

 void Anchors::init(AnchorType val, unsigned range)
{
    for (unsigned k = 0; k < range; k++)
        seqan::appendValue(set,val);
}

 void Anchors::init(int length)
{
    clear(set);
    seqan::appendValue(set, 0);
    (void)length;
}

 void Anchors::init()
{
    init(AnchorType(0), 0);
}

 void Anchors::setAnchor(unsigned p, 
    Anchors::AnchorType pos1,  Anchors::AnchorType pos2)
{
    set[p] = (pos1 << AnchorBase::bit) + pos2;
}

 Anchors::AnchorType Anchors::getPos1(unsigned p) const 
{
    return set[p] >> AnchorBase::bit;
}

 Anchors::AnchorType Anchors::getPos2(unsigned p) const
{
    return set[p] & AnchorBase::mask;
}

 Anchors::AnchorType Anchors::deltaPos1(unsigned p1, unsigned p2)
{
    return (set[p1] >> AnchorBase::bit) - (set[p2] >> AnchorBase::bit);
}

 Anchors::AnchorType Anchors::deltaPos2(unsigned p1, unsigned p2)
{
    return AnchorBase::mask & (set[p1] - set[p2]);
}
 Anchors::Iter Anchors::begin()
{
    return seqan::begin(set);
}
 Anchors::Iter Anchors::end()
{
    return seqan::end(set);
}
 void Anchors::sort(Anchors::Iter sortBegin, Anchors::Iter sortEnd)
{
    std::sort(sortBegin, sortEnd);
}
 void Anchors::sortPos2(Anchors::Iter sortBegin, Anchors::Iter sortEnd){
    AnchorBase::AnchorType mask = AnchorBase::mask;
    std::sort(sortBegin, sortEnd,
    [& mask](AnchorBase::AnchorType & a, 
                AnchorBase::AnchorType & b)
    {
        return (a & mask) < (b & mask);
    }) ;
}
 void Anchors::appendValue(Anchors::AnchorType val)
{
    seqan::appendValue(set,val);
}
uint64_t & Anchors::operator [](unsigned p)
{
    return set[p];
};
unsigned Anchors::length() 
{
    return seqan::length(set);
};

void MapParm::print()
{
    std::cerr << "blockSize " << blockSize << std::endl
            << "alpha " << alpha << std::endl
            << "alpha2 " << alpha2 << "\n"
            << "listN " << listN << "\n"
            << "listN2 " << listN2 << "\n"
            << "senThr " << senThr << "\n"
            << "delta " << delta << std::endl
            << "threshold " << threshold << std::endl
            << "kmerStep " << kmerStep << std::endl
            << "shapeLen " << shapeLen << std::endl
            //<<  "sensitivity " << sensitivity << "\n"
            << "anchorDeltaThr " << anchorDeltaThr << "\n"
            << "minReadLen " << minReadLen << "\n"
            << "anchorLenThr" << anchorLenThr << "\n"
            << "rcThr " << rcThr << "\n"
            << "cordThr" << cordThr << "\n";
}

static const String<Dna5> _complt = "tgcan";
 void _compltStr(String<Dna5> & str, String<Dna5> & res)
{
    resize(res, length(str));
    for (unsigned k = 0; k < length(str); k++)
     //   res[k]=_complt[str[k] - 'A'];
        res[k] = _complt[(unsigned)ordValue(str[k])];
}

 void _compltRvseStr(String<Dna5> & str, String<Dna5> & res)
{
    resize(res, length(str));
    for (unsigned k = 0; k < length(str); k++)
     //   res[k]=_complt[str[k] - 'A'];
    {
        res[k] = _complt[(unsigned)ordValue(str[length(str) - k - 1])];
        //std::cout << (unsigned)ordValue(str[length(str) - k - 1]) << std::endl;
    }
}

Dout dout;
Dout & Dout::operator << (int n)
{
    std::cout << n << " "; 
    return *this;
}
Dout & Dout::operator << (unsigned n)
{
    std::cout << n << " "; 
    return *this;
}
Dout & Dout::operator << (int64_t n)
{
    std::cout << n << " "; 
    return *this;
}
Dout & Dout::operator << (uint64_t n)
{
    std::cout << n << " "; 
}
Dout & Dout::operator << (CharString n)
{
    if (n == "\n")
    {
        std::cout << n;
    }
    else
    {
        std::cout << n << " "; 
    }
    return *this;
}
Dout & Dout::operator << (String<int64_t> & n)
{
    for (int i = 0; i < length(n); i++)
    {
        std::cout << n[i] << " ";
    }
    return * this;
}
Dout & Dout::operator << (double n)
{
    std::cout << n << " " ;
    return * this;
}

void ostreamWapper::print_message(std::string strs, 
                                  size_t start, 
                                  int end_type, 
                                  std::ostream & os)
{
    std::string spaces = "";
    if (start > length(contents))
    {
        for (int i = 0; i < start - length(contents); i++)
        {
            append (spaces, " ");
        }
    }
    append(contents, spaces);
    append(contents, strs);
    if (end_type == 1)
    {
        os << contents << "\n";
        contents = "";
    }
    else if (end_type == 2)
    {
        os << contents << "\r";
        contents = "";
    }
}

void ostreamWapper::print_message(double data, 
                                  size_t start, 
                                  int end_type, 
                                  std::ostream & os)
{
    float d = int(data * 100);
    std::ostringstream strs;
    strs << (d / 100);
    std::string str = strs.str();
    print_message(str, start, end_type, os);
}
void ostreamWapper::print_message(unsigned data, 
                                  size_t start, 
                                  int end_type, 
                                  std::ostream & os)
{
    std::ostringstream strs;
    strs << (data / 100);
    std::string str = strs.str();
    print_message(str, start, end_type, os);
}

ostreamWapper serr;

CmpInt64 & CmpInt64::init(int64_t & rslt, int64_t init_val) 
{ 
    p_rslt = & rslt;
    *p_rslt = init_val;
    return *this;
}
CmpInt64 & CmpInt64::min(int64_t & rslt, int64_t val) 
{
    return init(rslt, val);
}
CmpInt64 & CmpInt64::max(int64_t & rslt, int64_t val) 
{
    return init(rslt, val);
}
CmpInt64 & CmpInt64::operator << (int64_t n)
{
    if (*p_rslt > n)
    {
        *p_rslt = n;
    }
    return *this;
}
CmpInt64 & CmpInt64::operator >> (int64_t n)
{
    if (*p_rslt < n)
    {
        *p_rslt = n;
    }
    return *this;
}

