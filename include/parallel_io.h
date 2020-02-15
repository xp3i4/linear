#ifndef LINEAR_HEADER_PARALLEL_IO_H
#define LINEAR_HEADER_PARALLEL_IO_H

#include <seqan/sequence.h>
#include "f_io.h"

using seqan::String;

typedef String<Dna5> * P_Dna5s;
typedef CharString * P_CharStrings;
typedef String<uint64_t> * P_ULLs;
typedef String<BamAlignmentRecordLink> * P_BamLinks;

template<typename ElementType> //ElementType : P_Dna5s, etc.
class P_Buffer{
    String<ElementType> * p_buffer;
    int it_str;
    int it_end;
public:
    ElementType & operator [] (int i)
    {
        return (*p_buffer)[i];
    }
    int getBufferLength()
    {
        return length(*p_buffer);
    }
    int getEmptyBufferLength()
    {
        return it_end > it_str ? length(*p_buffer) - it_end + it_str : it_str - it_end;
    }
    void resize(int size)
    {
        seqan::resize (*p_buffer, size);   
    }
    void resize(int size, int val)
    {
        seqan::resize (*p_buffer, size, val);
    }
};

struct P_Tasks{
    /*
    P_Buffer<P_Dna5s> * p_p_reads;
    P_Buffer<P_CharStrings> * p_p_reads_ids;
    P_Buffer<P_ULLs> * p_p_cords; 
    P_Buffer<P_BamLinks> * p_p_bam_links;
    */
    std::list<StringSet<String<Dna5> > > * p_reads_buffer;
    std::list<StringSet<String<Dna5> > > * p_reads_ids_buffer;
    String<int> tasks;
    String<std::string> path1;
    String<std::string> path2;
    //critical resources
    int p_p_reads_str;
    int p_p_reads_end;
    int p_p_cords_str;
    int p_p_cords_end;
    int path1_i;
    int path2_i
    int sgn_request;
    int sgn_fetch;
    int sgn_print;
    SeqFileIn fin;
    std::ofstream fout;

    P_Tasks(P_Buffer<P_Dna5s> &, P_Buffer<P_CharStrings> &, P_Buffer<P_ULLs> &, P_Buffer<P_BamLinks> &, int);
    int setTasksFetchReads(int thread_id);
    int setTaskCalRecords(int thread_id);
    int setTasksPrintResults(int thread_id);
    int setTasksEnd(int thread_id);
    int isTasksFetchReads(int thread_id);
    int isTasksCalRecords(int thread_id);
    int isTasksPrintResults(int thread_id);
    int isTasksEnd(int thread_id);
    int isReadBufferFull();
    int isReadBufferEmpty();
};
struct  P_Parms
{

    int thd_read_block;

    P_Parms();
};

int p_FetchReads(P_Tasks & p_tasks);
int p_CalRecords();
int p_PrintResults();
int p_ThreadProcess(P_Tasks & p_tasks, P_Parms p_parms, int thread_id);

#endif
