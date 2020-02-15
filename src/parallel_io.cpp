#include "base.h"
#include "parallel_io.h"

using namespace seqan;

/*----------  Class P_Tasks  ----------*/
P_Tasks::P_Tasks(P_Buffer<P_Dna5s> & p_reads_buffer, P_Buffer<P_CharStrings> & p_reads_ids_buffer, P_Buffer<P_ULLs> & p_cords_buffer, P_Buffer<P_BamLinks> & p_bam_links_buffer, int thread_num
    String<std::string> g_paths, String<std::string> & r_paths)
{
    resize(tasks, thread_num, 1);
    p_p_reads = & p_reads_buffer;
    p_p_reads_ids = & p_reads_ids_buffer;
    p_p_cords = & p_cords_buffer;
    p_p_bam_links = & p_bam_links_buffer;
    paths1 = g_paths; 
    paths2 = r_paths;

    sgn_request = 0;
    sgn_fetch = 0;
    sgn_print = 0;
}
int P_Tasks::setTasksFetchReads(int thread_id)
{
    tasks[thread_id] = 2;
    return tasks[thread_id];
}
int P_Tasks::setTaskCalRecords(int thread_id)
{
    tasks[thread_id] = 4;
    return tasks[thread_id];
}
int P_Tasks::setTasksPrintResults(int thread_id)
{
    tasks[thread_id] = 8;
    return tasks[thread_id];
}
int P_Tasks::setTasksEnd(int thread_id)
{
    tasks[thread_id] == 0;
    return tasks[thread_id];
}
int P_Tasks::isTasksFetchReads(int thread_id)
{
    return tasks[thread_id] == 2;
}
int P_Tasks::isTasksCalRecords(int thread_id)
{
    return tasks[thread_id] == 4;
}
int P_Tasks::isTasksPrintResults(int thread_id)
{
    return tasks[thread_id] == 8;
}
int P_Tasks::isTasksEnd(int thread_id)
{
    return tasks[thread_id] == 0;
}
/*----------  Class P_Parms  ----------*/

P_Parms::P_Parms(String<std::string> & g_paths, String<std::string> & r_paths):
    thd_read_block(1000){}
/*----------  Threads functions  ----------*/

int p_RequestTask(P_Tasks & p_tasks, P_Parms & p_parms, int thread_id)
{
    //std::cerr << p_tasks.sgn_request << "\n";
    if (!atomicCas(p_tasks.sgn_request, 0, 1))
    {
            std::cerr << "xxx" << p_parms.path2[0] << "\n";
        P_Buffer<P_Dna5s> & pd = * p_tasks.p_p_reads;
        if (!atomicCas(p_tasks.sgn_fetch, 0, 1) &&  p_tasks.p_p_reads->getEmptyBufferLength() <= p_parms.thd_read_block)
        {
            p_tasks.setTasksFetchReads(thread_id);
        }
        /*
        else if (!atomicCas(p_tasks.sgn_print, 0, 1))
        {
            p_tasks.setTasksPrintResults(thread_id);
        }
        else if (!p_tasks.isReadBufferEmpty())
        {
            p_tasks.setTaskCalRecords(thread_id);
        }
        */
        atomicCas(p_tasks.sgn_request, 1, 0);
        return 1; 
    }
    else
    {
        return 0;
    }
}

int p_FetchReads(P_Tasks & p_tasks, P_Parms & p_parms)
{
    P_Buffer<P_Dna5s> & p_reads = *(p_tasks.p_p_reads);
    P_Buffer<P_CharStrings> & p_reads_ids = *(p_tasks.p_p_reads_ids);
    /*
    for (int i = p_parms.fetch_buffer_str; i < p_parms.fetch_buffer_end; i++)
    {
        P_Dna5s p_read = new String<Dna5>;
        P_CharStrings p_read_id = new CharString;
        readRecord (*p_read, *p_read_id, file);
        p_reads[i] = p_read;
        p_reads_ids[i] = p_read_id;
    }
    */
    return 0;
}
/*
 * Calculate records
 */
int p_CalRecords()
{

    return 0;
}

int p_PrintResults()
{
    return 0; 
}

int p_ThreadProcess(P_Tasks & p_tasks, P_Parms p_parms, int thread_id)
{
    while (!p_RequestTask(p_tasks, p_parms, thread_id) && !p_tasks.isTasksEnd(thread_id))
    {
        if(p_tasks.isTasksFetchReads(thread_id)) 
        {
            p_FetchReads(p_tasks, p_parms);
        }
        /*
        else if (p_tasks.isTasksCalRecords(thread_id))
        {
            p_CalRecords();
        }
        else if (p_tasks.isTasksPrintResults(thread_id))
        {
            p_PrintResults();
        }
        */
    } 
    return 0;
}
