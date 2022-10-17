#include "base.h"
#include "ska_sort.hpp"
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

unsigned const UMAX = (1 << 30) - 1;
uint64_t const ULLMAX = ~0;
int64_t const LLMAX = (1LL << 62) - 1; //(1ULL << 63) - 1 integer overflow on some compilers
int64_t const LLMIN = -LLMAX;

bool f_debug = true;

using std::cerr;

Options::Options():
        op_status(0),
        name("Linear"),
        version("1.8.2"),
        //slogan("\033[1;31mN\033[movel \033[1;34mE\033[mfficient \033[1;33mC\033[moncise"),
        //slogan("A\033[1;31mL\033[m\033[1;32mI\033[mg\033[1;33mN\033[mment-fre\033[1;34mE\033[m method for long-read v\033[1;35mA\033[m\033[1;36mR\033[miants resolution"),
        slogan("Extensible Long-read Algorithms Framework"),
        //oPath(""),
        gap_len(1),
        apx_chain_flag(1),
        aln_flag(0),
        sam_flag(1),
        apf_flag(1),
        reform_ccs_cigar_flag(0),
        bal_flag(1),
        f_output_type(2),
        f_dup(0),
        sensitivity(1),
        thread(16),
        index_t(1),
        feature_t(2), //apx2_48 by defalut
        read_group(""),
        sample_name(""),
        sequence_sam_flag(1)
        {
           date += __TIME__; 
           date += " ";
           date += __DATE__;
        }
std::string Options::getOutputPath() const {return oPath;}
void Options::printRunInfo(){
    std::cerr << name << ": " << slogan << std::endl; 
}

Options::Options(int argc, char const ** & argv) : Options()
{
    op_argc = argc;
    op_argv = argv;
    if (length(argv) < 1)
    {
        append(cmd_line, CharString(argv[1]));
        for (int i = 1; i < argc; i++)
        {
            append(cmd_line, " ");
            append(cmd_line, CharString(argv[i]));
        }
    }
}


unsigned  Options::isOutputApf()
{
    return f_output_type & 1;
}
unsigned Options::isOutputSam()
{
    return f_output_type & 2;
}

unsigned Options::isOutputBamPbsv()
{
    return f_output_type & 8;
}
unsigned Options::isOutputBamStandard()
{
    return f_output_type & 4;
}

void Parms::toggle(int i)
{(void) i;}

GlobalParms::GlobalParms()
{
    shape_len = 21;
}


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

std::ifstream::pos_type _filesize(const char* filename)
{
    std::ifstream in(filename,  std::ifstream::binary | std::ifstream::ate);
    return in.tellg(); 
}

/*
 *[]::load all genome from path [.fa]
 */
std::pair<uint, uint> loadRecords(StringSet<String<Dna5> > & seqs, 
            StringSet<CharString> & ids, 
            Options::PathType path)
{
    //double time = sysTime();
    SeqFileIn gFile;
    std::pair<uint, uint> res;
    if (!open(gFile, toCString(path)))
    {
        serr.print_message("\033[1;31mE[03]:\033[0m can't open file ", 0, 0, std::cerr);
        serr.print_message(toCString(path), 0, 1, std::cerr);
        res =std::make_pair (uint(~0), uint(~0));
        return res;
    }
    //double fileSize = _filesize (toCString(path));
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
                unsigned k = 1;
                while (!flag)
                {
                    serr.print_message("", 50, 2, std::cerr);
                    std::cerr << "=>Read genomes " << seqCount << dotstatus[(k - 1)/10 %3] << "\r"; 
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
                        serr.print_message("\033[1;31mE[04]:\033[0m can't read geonmes in file", 0, 1, std::cerr);
                        error = 1;
                    }
                    for (unsigned i = 0; i < length(tmp_id); i++)
                    {
                        if (tmp_id[i] == ' ')
                        {
                            resize(tmp_id, i);
                            break;
                        }
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
        if (len_sum == uint(~0) || seqCount == uint(~0))
        {
            status += 1;
            continue;
        }
        if (i == 0)
        {
            serr.print_message ("--Read genomes ", 0, 0, cerr);
            serr.print_message (" ", 50, 1, cerr);
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
    (void) options;
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
    //std::sort(sortBegin, sortEnd);
    ska_sort(sortBegin, sortEnd);
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
}
unsigned Anchors::length() 
{
    return seqan::length(set);
}

static const String<Dna5> _complt = "tgcan";
 void _compltStr(String<Dna5> & str, String<Dna5> & res)
{
    resize(res, length(str));
    for (unsigned k = 0; k < length(str); k++)
    {
        res[k] = _complt[(unsigned)ordValue(str[k])];
    }
}

 void _compltRvseStr(String<Dna5> & str, String<Dna5> & res)
{
    resize(res, length(str));
    for (unsigned k = 0; k < length(str); k++)
    //res[k]=_complt[str[k] - 'A'];
    {
        res[k] = _complt[(unsigned)ordValue(str[length(str) - k - 1])];
        //std::cout << (unsigned)ordValue(str[length(str) - k - 1]) << std::endl;
    }
}
/*----------  Debug ostream to replace ostream ----------*/
Gout dout;
Gout::Gout(bool f_p)
{
    f_print = f_p;
}
/*
Gout & operator << (Dout & d, int n)
{
    unused(d);
    Gout *p = new Gout;
    *p << n;
    return *p;
}
Gout & operator << (Dout & d, unsigned n)
{
    unused(d);
    Gout *p = new Gout;
    *p << n;
    return *p;  
}
Gout & operator << (Dout & d, int64_t n)
{
    unused(d);
    Gout *p = new Gout;
    *p << n;
    return *p;
}
Gout & operator << (Dout & d, uint64_t n)
{
    unused(d);
    Gout *p = new Gout;
    *p << n;
    return *p;
}
Gout & operator << (Dout & d, CharString n)
{
    unused(d);
    Gout *p = new Gout;
    *p << n;
    return *p;
}
Gout & operator << (Dout & d, String<int64_t> & n)
{
    unused(d);
    Gout *p = new Gout;
    *p << n;
    return *p;
}
Gout & operator << (Dout & d, double n)
{
    unused(d);
    Gout *p = new Gout;
    *p << n;
    return *p;
}
*/
Gout & Gout::operator << (int n)
{
    if (f_print)
    {
        buffer << n << " "; 
    }
    return *this;
}
Gout & Gout::operator << (unsigned n)
{
    if (f_print)
    {
        buffer << n << " "; 
    }
    return *this;
}
Gout & Gout::operator << (int64_t n)
{
    if (f_print)
    {
        buffer << n << " "; 
    }
    return *this;
}
Gout & Gout::operator << (uint64_t n)
{
    if (f_print)
    {
        buffer << n << " "; 
    }
    return *this;
}
Gout & Gout::operator << (CharString n)
{
    if (f_print)
    {
        buffer << n << " ";
        if (n == "\n")
        {
            std::cout << buffer.str();
            buffer.str("");
        }
    }
    return *this;
}
Gout & Gout::operator << (String<int64_t> & ns)
{
    if (f_print)
    {
        for (auto n : ns)
        {
            buffer << n << " ";
        }
    }  
    return * this;
}
Gout & Gout::operator << (double n)
{
    if (f_print)
    {
        buffer << n << " " ;
    }
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
        for (uint i = 0; i < start - length(contents); i++)
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
    strs << data;
    std::string str = strs.str();
    print_message(str, start, end_type, os);
}
void ostreamWapper::print_message(int data, 
                                  size_t start, 
                                  int end_type, 
                                  std::ostream & os)
{
    std::ostringstream strs;
    strs << data;
    std::string str = strs.str();
    print_message(str, start, end_type, os);
}

ostreamWapper serr;

std::string _2stdstring (CharString str)
{
    std::string rsl;
    for (unsigned i = 0; i < length(str); i++)
    {
        rsl.push_back(char(str[i]));
    }
    return rsl;
}

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

void sort_ska(Iterator<String<uint64_t> >::Type it_str, Iterator<String<uint64_t> >::Type it_end)
{
    ska_sort(it_str, it_end);
}

int print_seq(String<Dna5> & seq, uint64_t str, uint64_t end, std::string header)
{
    std::cout << header << " ";
    for (unsigned i = str; i < end && i < length(seq); i++)
    {
        std::cout << seq[i];
    }
    std::cout << "\n";
    return 0;
}
int mod(int a, int b){int c = a % b; return c >= 0 ? c : c + b;}



void F_Print_::setPrintApf(uint & f){f |= 1;}
void F_Print_::unsetPrintApf(uint & f){f &= ~1;}
void F_Print_::setPrintSam(uint & f){f |= 2;}
void F_Print_::unsetPrintSam(uint & f){f &= ~2;}
void F_Print_::setPrintBamStd(uint & f){f |= 4;}
void F_Print_::unsetPrintBamStd(uint & f){f &= ~4;}
void F_Print_::setPrintBamPbsv(uint & f){f |= 8;}
void F_Print_::unsetPrintBamPbsv(uint & f){f &= ~8;}
void F_Print_::clear(uint & f){f = 0;}
int F_Print_::isPrintApf(uint f){return f & 1;}
int F_Print_::isPrintSam(uint f){return f & 2;}
int F_Print_::isPrintBam(uint f){return isPrintBamPbsv(f) || isPrintBamStd(f);}
int F_Print_::isPrintBamStd(uint f){return f & 4;}
int F_Print_::isPrintBamPbsv(uint f){return f & 8;}
int F_Print_::isPrintSamBam(uint f){return isPrintSam(f) || isPrintBam(f);}

F_Print_ fp_handler_;
