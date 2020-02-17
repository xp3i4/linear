#include <iostream>
#include <seqan/stream.h>
#include "base.h"
#include "parallel_io.h"

using namespace seqan;
/*----------  Class P_Task  ----------*/
P_Task::P_Task()
{
    //0:task end; 1:normal(wait); 2:fetch reads; 4:calculate Records; 8:print results 
    task_type = 1;
}
/*----------  Class P_Tasks  ----------*/
P_Tasks::P_Tasks(P_ReadsBuffer & reads_buffer, P_ReadsIdsBuffer & reads_ids_buffer, 
    StringSet<std::string> & g_paths, StringSet<std::string> & r_paths, int thread_num)
{
    resize(tasks, thread_num);
    for (int i = 0; i < length(tasks); i++)
    {
        //tasks[i].task_type = 0;
    }
    p_reads_buffer = & reads_buffer;
    p_reads_ids_buffer = & reads_ids_buffer;

    paths1 = g_paths; 
    paths2 = r_paths;

    paths1_it = 0;
    paths2_it = 0;
    sgn_request = 0;
    sgn_fetch = 0;
    sgn_fetch_end = 0;
    sgn_print = 0;
    sgn_all_tasks_end = 0;
    f_fin_open = 0;
}
void P_Tasks::setTasksWait(int thread_id)
{
    tasks[thread_id].task_type = 1;
}
void P_Tasks::setTasksFetchReads(int thread_id)
{
    tasks[thread_id].task_type = 2;
}
void P_Tasks::setTaskCalRecords(int thread_id)
{
    tasks[thread_id].task_type = 4;
}
void P_Tasks::setTasksPrintResults(int thread_id)
{
    tasks[thread_id].task_type = 8;
}
void P_Tasks::setTasksEnd(int thread_id)
{
    tasks[thread_id].task_type = 0;
}
void P_Tasks::setTasksAllEnd()
{
    atomicCas(sgn_all_tasks_end, 0, 1);
}
int P_Tasks::isTasksFetchReads(int thread_id)
{
    return tasks[thread_id].task_type == 2;
}
int P_Tasks::isTasksCalRecords(int thread_id)
{
    return tasks[thread_id].task_type == 4;
}
int P_Tasks::isTasksPrintResults(int thread_id)
{
    return tasks[thread_id].task_type == 8;
}
int P_Tasks::isTasksEnd(int thread_id)
{
    return tasks[thread_id].task_type == 0;
}
int P_Tasks::isTasksAllEnd()
{
    return atomicCas(sgn_all_tasks_end, 1, sgn_all_tasks_end); //just get the value atomically 
}
/*----------  Class P_Parms  ----------*/

P_Parms::P_Parms():
    thd_fetch_num(1),
    thd_fetch_block_size(5)
    {}
/*----------  Threads functions  ----------*/
/*
 * Free critical variables 
 */
int freeCurrentTask_(P_Tasks & p_tasks, P_Parms & p_parms, int thread_id, int f_error)
{
    if (f_error)
    {
        p_tasks.setTasksAllEnd();
    }
    if (p_tasks.isTasksFetchReads(thread_id))
    {
        atomicCas(p_tasks.sgn_fetch, 1, 0);
        p_tasks.setTasksWait(thread_id);
    }
    else if (p_tasks.isTasksCalRecords(thread_id))
    {

    }
    else if (p_tasks.isTasksPrintResults(thread_id))
    {

    }
    else
    {
        //none
    }
    return 0;
}
/*
 * Called after calling function @freeCurrentTask_
 */
int requestNewTask_(P_Tasks & p_tasks, P_Parms & p_parms, int thread_id, int f_error)
{
    if (p_tasks.isTasksAllEnd())  
    {
        p_tasks.setTasksEnd(thread_id);
    }
    else if (p_tasks.sgn_fetch_end != 1)
    {
        if (!atomicCas(p_tasks.sgn_fetch, 0, 1))
        {
            p_tasks.setTasksFetchReads(thread_id);
        }
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
    else
    {
        if (p_tasks.sgn_fetch_end == 1) // !!!Add all end conditions here !!!
        {
            p_tasks.setTasksEnd(thread_id);
        }
    }   
    return 0;
}
/*
 * Tasks main controller
 */
int p_RequestTask(P_Tasks & p_tasks, P_Parms & p_parms, int thread_id, int f_error)
{
    //std::cerr << p_tasks.sgn_request << "\n";
    if (!atomicCas(p_tasks.sgn_request, 0, 1))
    {
        freeCurrentTask_(p_tasks, p_parms, thread_id, f_error);
        requestNewTask_ (p_tasks, p_parms, thread_id, f_error);
        atomicCas(p_tasks.sgn_request, 1, 0);
        return 0; 
    }
    return 1;
}

int p_FetchReads(P_Tasks & p_tasks, P_Parms & p_parms)
{
    P_ReadsIdsBuffer & buffer1 = *p_tasks.p_reads_ids_buffer;
    P_ReadsBuffer & buffer2 = *p_tasks.p_reads_buffer;
    if (p_tasks.paths2_it >= length(p_tasks.paths2))
    {
        p_tasks.sgn_fetch_end = 1;
        return 1;
    }
    std::string file_name = p_tasks.paths2[p_tasks.paths2_it];
    for (int i = 0; i < p_parms.thd_fetch_num && ! buffer1.isFull(); i++)
    {
        if (!p_tasks.f_fin_open)
        {
            if(!open(p_tasks.fin, toCString(file_name)))
            {
                serr.print_message("\033[1;31mError:\033[0m can't open read file ", 2, 0, std::cerr);
                serr.print_message(toCString(file_name), 0, 1, std::cerr);
                return 2; 
            }  
            else
            {
                p_tasks.f_fin_open = 1;
            }  
        }
        clear(buffer1.in());
        clear(buffer2.in());
        readRecords(buffer1.in(), buffer2.in(), p_tasks.fin, p_parms.thd_fetch_block_size);
        //<<debug
        /*
        for (int i = 0; i < length(buffer1.out()); i++)
        {
            dout << "p_se " << i <<  length(buffer2.out()[i]) << buffer2.usedSize() << length(buffer1.out()) << "\n";
        }
        buffer1.nextOut();
        buffer2.nextOut();
        */
        //>>debug
        buffer1.nextIn();
        buffer2.nextIn();
    }
    if (atEnd(p_tasks.fin))
    {
        close(p_tasks.fin);
        p_tasks.paths2_it++;
        p_tasks.f_fin_open = 0;
    }
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
        //std::cout << "p1" << "\n";
    int f_error = 0;
    while (true)
    {
        if (p_RequestTask(p_tasks, p_parms, thread_id, f_error))
        {   
            continue; //continue request if failed in requesting task
        }
        if (p_tasks.isTasksEnd(thread_id))
        {
            break;
        }
        else if(p_tasks.isTasksFetchReads(thread_id)) 
        {
            f_error = p_FetchReads(p_tasks, p_parms);
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
