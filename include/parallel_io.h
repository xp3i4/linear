#ifndef LINEAR_HEADER_PARALLEL_IO_H
#define LINEAR_HEADER_PARALLEL_IO_H

#include <seqan/sequence.h>
#include "f_io.h"
/*
 * This is a simple model for parallel I/O buffers.
 * Restrictions : Only one thread is allowed to input or output at the same time,
   while parallel of I/O and jobs are allowed. 
 * Namely one thread fetch read from the disk, 
   multiple threads for jobs and one thread dump results to the disk simultaneously!
 */
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
//NOTE::!!
//Functions of P_Buffer are not thread safe for reasons regarding computational performance.
//Use them carefully!
template<class T>
class P_Buffer
{
    //use size() + 1 length(physical len) of list to simulate size() length (virtual len) of buffer
    //to avoid ambiguity of both empty and full has it_in == it_out if phy len == virtual len 
    RandomList<T> buffer;
    String<int> buffer_status;
    int it_in;
    int it_out;
    int CONST_F_PROTECTED; //it_in or it_out is not allowed to access/skip protected variables;
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
typedef P_Buffer<std::string> P_ReadsPathsBuffer;
typedef P_Buffer<StringSet<String<uint64_t> > > P_CordsBuffer;
typedef P_Buffer<StringSet<String<CordInfo> > > P_CordsInfoBuffer;
typedef P_Buffer<StringSet<String<BamAlignmentRecordLink> > > P_BamLinkBuffer;

struct P_Parms
{
    int thd_fetch_num;
    int thd_assign_num;
    int thd_print_num;
    int thd_buffer_block_size;
    P_Parms(int, int, int, int);
    P_Parms();
    void printParms();
};

struct Counters{
    String<int> counters;
    String<double> timers;
    //add here the new counters 
    Counters(){resize(counters, 3, 0); resize(timers, 3, 0);}
    int getInCounter(){return counters[0];}
    int getCalCounter(){return counters[1];}
    int getOutCounter(){return counters[2];}
    double getInTimer(){return timers[0];}
    double getCalTimer(){return timers[1];}
    double getOutTimer(){return timers[2];}
    void setInCounter(int val){counters[0] = val;}
    void setCalCounter(int val){counters[1] = val;}
    void setOutCounter(int val){counters[2] = val;}
    void setInTimer(double val){timers[0] = val;}
    void setCalTimer(double val){timers[1] = val;}
    void setOutTimer(double val){timers[2] = val;}
    
    void clearCounters(){
        counters[0] = counters[1] = counters[2] = 0;
        timers[0] = timers[1] = timers[2] = 0;
    }
    void addCounters (Counters & b_counters)
    {
        counters[0] += b_counters.counters[0];
        counters[1] += b_counters.counters[1];
        counters[2] += b_counters.counters[2];
        timers[0] += b_counters.timers[0];
        timers[1] += b_counters.timers[1];
        timers[2] += b_counters.timers[2];
    }
};

struct P_Task{
    //all input buffers are supposed to have the same length. This applies to the output buffers too.
    String<int> p_ins; //points to i_th element in buffer(input) for calculation
    String<int> p_outs; //points to cords and bam for buffering result during the calculation
    Counters counters;
    int task_type;
    int thread_id;
    P_Task();
};

struct P_Tasks{
    double start__time;
    //buffers
    //Warn!!::make sure the in and out iterator of the Group are modified 
    //by at most one function at the same time
    //Group1:buffers of Input
    P_ReadsBuffer * p_reads_buffer;
    P_ReadsIdsBuffer * p_reads_ids_buffer;
    P_ReadsPathsBuffer * p_reads_paths_buffer;
    //Group2:buffers of output
    P_CordsBuffer * p_cords1_buffer; //cords_str
    P_CordsBuffer * p_cords2_buffer; //cords_end
    P_CordsInfoBuffer * p_cords_info_buffer;
    P_BamLinkBuffer * p_bam_link_buffer;

    //tasks of each thread : tasks[thread_id] 
    String<P_Task> tasks;
    Counters counters;

    //critical resources
    StringSet<std::string> paths1;
    StringSet<std::string> paths2;
    int paths1_it;
    int paths2_it;
    /*WARN::assign_it only be used in P_Tasks::assignCalRecords; To prevent 
      synchronization conflicat, Don't use them in other functions.*/
    int assign_it2; //iter to read buffer ready for assignment
    int assign_it3; //iter to cords1 2 and bam_buffer
    int sgn_request;
    int sgn_fetch;
    int sgn_fetch_end;
    int sgn_print;
    int sgn_all_tasks_end;
    int f_fin_open;
    int f_printRunningInfos;
    seqan::SeqFileIn fin;
    std::ofstream fout;
    void initPTasks(P_ReadsBuffer &,  P_ReadsIdsBuffer &, P_ReadsPathsBuffer &,
                    P_CordsBuffer &, P_CordsBuffer &, P_CordsInfoBuffer &, P_BamLinkBuffer &,
        StringSet<std::string> &, StringSet<std::string> &, int);
    P_Tasks();
    P_Tasks(P_ReadsBuffer &, P_ReadsIdsBuffer &, P_ReadsPathsBuffer&, 
        P_CordsBuffer &, P_CordsBuffer &, P_CordsInfoBuffer &, P_BamLinkBuffer &, 
        StringSet<std::string> &, StringSet<std::string> &, int);
    void startRunTasks();
    //printInfos
    void printRunningInfos();
    //iterator
    void nextAssignIt2(){p_reads_buffer->nextIt(assign_it2);}
    void nextAssignIt3(){p_cords1_buffer->nextIt(assign_it3);}
    void nextGroup1In();
    void nextGroup2In();
    //modifier
    void setTaskWait (int thread_id);
    void setTaskFetchReads(int thread_id);
    void setTaskCalRecords(int thread_id);
    void setTaskPrintResults(int thread_id);
    void setTaskEnd(int thread_id);
    void setTaskAllEnd();
    void updateCounters();
    //is
    int isTaskFetchReads(int thread_id);
    int isTaskCalRecords(int thread_id);
    int isTaskPrintResults(int thread_id);
    int isTaskEnd(int thread_id);
    int isTaskAllEnd();
    int isAllInBuffsEmpty();
    int isAllOutBuffsEmpty();
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

//Wrapper class to abstract two operations
class P_Mapper
{
protected:
    P_Tasks p_tasks;
    P_Parms p_parms;
    P_ReadsBuffer reads_buffer;
    P_ReadsIdsBuffer reads_ids_buffer; 
    P_ReadsPathsBuffer reads_paths_buffer;
    P_CordsBuffer cords1_buffer;
    P_CordsBuffer cords2_buffer;
    P_CordsInfoBuffer cords_info_buffer;
    P_BamLinkBuffer bam_link_buffer;

public:
    P_Task & getPTask(int);
    P_Tasks & getPTasks();
    P_ReadsBuffer & getPReadsBuffer();
    P_ReadsIdsBuffer & getPReadsIdBuffer();
    P_ReadsPathsBuffer & getPReadsPathsBuffer();
    P_CordsBuffer & getPCords1Buffer();
    P_CordsBuffer & getPCords2Buffer();
    P_CordsInfoBuffer & getPCordsInfoBuffer();
    P_BamLinkBuffer & getPBamLinksBuffer();
    void initBuffers(int reads_buffer_size, int cords_buffer_size, 
        StringSet<std::string> &, StringSet<std::string> &, int, P_Parms &);

    virtual int p_calRecords(int, int, int){return 0;}; //input: in_buffer_id, out_buffer_id
    virtual int p_printResults(int, int, int){return 0;};
};

int p_ThreadProcess(P_Mapper & mapper, P_Parms & p_parms, int thread_id);
/*----------  Global utilities  ----------*/
//None

#endif
