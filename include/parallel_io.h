#ifndef LINEAR_HEADER_PARALLEL_IO_H
#define LINEAR_HEADER_PARALLEL_IO_H

#include <seqan/sequence.h>
#include "f_io.h"

using seqan::String;
/*----------  A wrapper of Random access list  ----------*/

template<class T>
class RandomList{
    typedef typename std::list<T>::iterator RLIter;
    std::list<T> r_list;
    String<RLIter> its;
public:
    T & operator [] (int i)
    {
        return *its[i];
    }
    int size(){
        return r_list.size();
    }
    void resize(int len)
    {
        r_list.resize(len);
        seqan::resize(its, len);
        RLIter it = r_list.begin();
        for (int i = 0; i < len; i++)
        {
            its[i] = it++;
        }
    }
    RandomList(){}
};
/*----------  A wrapper of buffer  ----------*/

template<class T>
class P_Buffer
{
    //use size() + 1 length(physical len) of list to simulate size() length (virtual len) of buffer
    //to avoid ambiguity of both empty and full has it_in == it_out if phy len == virtual len 
    RandomList<T> buffer;
    String<int> buffer_status;
    int it_in;
    int it_out;
    int CONST_F_PROTECTED; //it_in or it_out can't access or skip protected variables;
public:
    T & operator [] (int i){return buffer[i];}
    T & in() {return buffer[it_in];}    //return the T for input
    T & out() {return buffer[it_out];}  //..output
    int inIt(){return it_in;}
    int outIt(){return it_out;}
    int physicalSize(){return buffer.size();}
    int size() {return this->physicalSize() - 1;}  //return virtual size 
    int size(int it_str, int it_end){
        return mod((it_end - it_str), physicalSize());
    } //size of [it_str, it_end)
    int usedSize(){return mod((it_in - it_out), this->physicalSize());}
    int isFull() {return usedSize() == this->size();}
    int isEmpty() {return usedSize() == 0;}
    int isProtected(int i){return buffer_status[i] & CONST_F_PROTECTED;}
    //modifer:=== 
    void setProtected(int i){buffer_status[i] |= CONST_F_PROTECTED;}
    void unsetProtected(int i){buffer_status[i] &= (~ CONST_F_PROTECTED);}
    void resize(int len) 
    {
        buffer.resize(len + 1);
        seqan::resize(buffer_status, len + 1, CONST_F_PROTECTED);
    } // len is the virtual size
    void nextIt(int & it){it = mod(++it, this->physicalSize());}
    void nextOut(){nextIt(it_out);}
    void nextIn(){nextIt(it_in);} 
    P_Buffer()
    {
        it_in = 0; it_out = 0; 
        CONST_F_PROTECTED = 1;
    }
};

/*----------  Parallel tasks controller  ----------*/

typedef StringSet<String<Dna5> > * P_P_Reads;
typedef P_Buffer<StringSet<String<Dna5> > > P_ReadsBuffer;
typedef P_Buffer<StringSet<CharString> > P_ReadsIdsBuffer;
typedef P_Buffer<StringSet<String<uint64_t> > > P_CordsBuffer;
typedef P_Buffer<StringSet<String<BamAlignmentRecordLink> > > P_BamLinkBuffer;

struct P_Parms
{
    int thd_fetch_num;
    int thd_buffer_block_size;
    int thd_assign_num;
    int thd_print_num;
    P_Parms();
    void printParms();
};

struct P_Task{
    //all input buffers are supposed to have the same length. This applies to the output buffers too.
    String<int> p_ins; //points to i_th element in buffer(input) for calculation
    String<int> p_outs; //points to cords and bam for output result
    int task_type;
    int thread_id;
    P_Task();
};

struct P_Tasks{
    //buffers
    P_ReadsBuffer * p_reads_buffer;
    P_ReadsIdsBuffer * p_reads_ids_buffer;
    P_CordsBuffer * p_cords_buffer;
    P_BamLinkBuffer * p_bam_link_buffer;

    //settings of task of each thread
    String<P_Task> tasks;

    //critical resources
    StringSet<std::string> paths1;
    StringSet<std::string> paths2;
    int paths1_it;
    int paths2_it;
    int assign_it2; //iter to read ready for assignment for calculatation
    int assign_it3;
    int sgn_request;
    int sgn_fetch;
    int sgn_fetch_end;
    int sgn_print;
    int sgn_all_tasks_end;
    int f_fin_open;
    seqan::SeqFileIn fin;
    std::ofstream fout;

    P_Tasks(P_ReadsBuffer &, P_ReadsIdsBuffer &, P_CordsBuffer &, P_BamLinkBuffer &, StringSet<std::string> &, StringSet<std::string> &, int);
    //modifier
    void setTaskWait (int thread_id);
    void setTaskFetchReads(int thread_id);
    void setTaskCalRecords(int thread_id);
    void setTaskPrintResults(int thread_id);
    void setTaskEnd(int thread_id);
    void setTaskAllEnd();
    //is
    int isTaskFetchReads(int thread_id);
    int isTaskCalRecords(int thread_id);
    int isTaskPrintResults(int thread_id);
    int isTaskEnd(int thread_id);
    int isTaskAllEnd();
    //int isReadBufferFull();
    //int isReadBufferEmpty();
    //
    int getTaskInBufferLen(int thread_id);
    int getTaskOutBufferLen(int thread_id);
    int getTaskInBufferPtr(int thread_id, int i);
    int getTaskOutBufferPtr(int thread_id, int i);
    //
    int assignCalRecords(P_Parms & p_parms, int thread_id);
};


int p_FetchReads(P_Tasks & p_tasks);
int p_CalRecords();
int p_PrintResults();
int p_ThreadProcess(P_Tasks & p_tasks, P_Parms p_parms,

                    int thread_id);
//Wrapper class to abstract two operations
class P_Mapper
{
protected:
    P_ReadsBuffer reads_buffer;
    P_ReadsIdsBuffer reads_ids_buffer; 
    P_CordsBuffer cords_buffer;
    P_BamLinkBuffer bam_link_buffer;
public:
    calRecords();
    printResults();
}

/*----------  Global utilities  ----------*/
//None

#endif
