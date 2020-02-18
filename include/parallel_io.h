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
    int it_in;
    int it_out;
public:
    T & operator [] (int i){return buffer[i];}
    T & in() {return buffer[it_in];}    //return the T for input
    T & out() {return buffer[it_out];}  //..output
    int InIt(){return it_in;}
    int OutIt(){return it_out;}
    int physicalSize(){return buffer.size();}
    int size() {return this->physicalSize() - 1;}  //return virtual size 
    int size(int it_str, int it_end){
        return (it_end - it_str) % physicalSize();
    } //size of [it_str, it_end)
    int usedSize(){return (it_in - it_out) % this->physicalSize();}
    int isFull() {return usedSize() == this->size();}
    int isEmpty() {return usedSize() == 0;}
    //modifer: 
    void resize(int len) {buffer.resize(len + 1);} // len is the virtual size
    void nextOut(){it_out = ++it_out % this->physicalSize();}
    void nextIn(){int p = it_in; it_in = ++it_in % this->physicalSize();dout << "nextin---" << p << it_in << "\n";} 
    P_Buffer(){it_in = 0; it_out = 0;}
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
    int thd_fetch_block_size;
    int thd_assign_num;
    P_Parms();
};

struct P_Task{
    String<int> p_reads; //*p_p_reads[i] to get the StringSet<..>
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
    int assign_it; //iter to read ready for assignment for calculatation
    int sgn_request;
    int sgn_fetch;
    int sgn_fetch_end;
    int sgn_print;
    int sgn_all_tasks_end;
    int f_fin_open;
    seqan::SeqFileIn fin;
    std::ofstream fout;

    P_Tasks(P_ReadsBuffer &, P_ReadsIdsBuffer &, StringSet<std::string> &, StringSet<std::string> &, int);
    void setTasksWait (int thread_id);
    void setTasksFetchReads(int thread_id);
    void setTaskCalRecords(int thread_id);
    void setTasksPrintResults(int thread_id);
    void setTasksEnd(int thread_id);
    void setTasksAllEnd();
    int isTasksFetchReads(int thread_id);
    int isTasksCalRecords(int thread_id);
    int isTasksPrintResults(int thread_id);
    int isTasksEnd(int thread_id);
    int isTasksAllEnd();
    int assignCalReads(P_Parms & p_parms, int thread_id);
    //int isReadBufferFull();
    //int isReadBufferEmpty();
};


int p_FetchReads(P_Tasks & p_tasks);
int p_CalRecords();
int p_PrintResults();
int p_ThreadProcess(P_Tasks & p_tasks, P_Parms p_parms, int thread_id);

/*----------  Global utilities  ----------*/
//None

#endif
