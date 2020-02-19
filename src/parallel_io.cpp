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
    P_CordsBuffer & cords_buffer, P_BamLinkBuffer & bam_link_buffer, 
    StringSet<std::string> & g_paths, StringSet<std::string> & r_paths, int thread_num)
{
    resize(tasks, thread_num);
    for (int i = 0; i < length(tasks); i++)
    {
        //tasks[i].task_type = 0;
    }
    p_reads_buffer = & reads_buffer;
    p_reads_ids_buffer = & reads_ids_buffer;
    p_cords_buffer = & cords_buffer;
    p_bam_link_buffer = & bam_link_buffer;

    paths1 = g_paths; 
    paths2 = r_paths;

    paths1_it = 0;
    paths2_it = 0;
    assign_it2 = p_reads_buffer->outIt(); // = 0;
    assign_it3 = p_cords_buffer->outIt(); // = 0;
    sgn_request = 0;
    sgn_fetch = 0;
    sgn_fetch_end = 0;
    sgn_print = 0;
    sgn_all_tasks_end = 0;
    f_fin_open = 0;
}
void P_Tasks::setTaskWait(int thread_id)
{
    tasks[thread_id].task_type = 1;
}
void P_Tasks::setTaskFetchReads(int thread_id)
{
    tasks[thread_id].task_type = 2;
}
void P_Tasks::setTaskCalRecords(int thread_id)
{
    tasks[thread_id].task_type = 4;
}
void P_Tasks::setTaskPrintResults(int thread_id)
{
    tasks[thread_id].task_type = 8;
}
void P_Tasks::setTaskEnd(int thread_id)
{
    tasks[thread_id].task_type = 0;
}
void P_Tasks::setTaskAllEnd()
{
    atomicCas(sgn_all_tasks_end, 0, 1);
}
int P_Tasks::isTaskFetchReads(int thread_id)
{
    return tasks[thread_id].task_type == 2;
}
int P_Tasks::isTaskCalRecords(int thread_id)
{
    return tasks[thread_id].task_type == 4;
}
int P_Tasks::isTaskPrintResults(int thread_id)
{
    return tasks[thread_id].task_type == 8;
}
int P_Tasks::isTaskEnd(int thread_id)
{
    return tasks[thread_id].task_type == 0;
}
int P_Tasks::isTaskAllEnd()
{
    return atomicCas(sgn_all_tasks_end, 1, sgn_all_tasks_end); //just get the value atomically 
}
int P_Tasks::getTaskInBufferPtr(int thread_id, int i)
{
    return tasks[thread_id].p_ins[i];
}
int P_Tasks::getTaskOutBufferPtr(int thread_id, int i)
{
    return tasks[thread_id].p_outs[i];
}
int P_Tasks::getTaskInBufferLen(int thread_id)
{
    return length(tasks[thread_id].p_ins);
}
int P_Tasks::getTaskOutBufferLen(int thread_id)
{
    return length(tasks[thread_id].p_outs);
}
int P_Tasks::assignCalRecords(P_Parms & p_parms, int thread_id)
{
    if (!empty(tasks[thread_id].p_ins) || !empty(tasks[thread_id].p_outs))
    {
        return 1;
    }
    P_ReadsBuffer &     buffer2 = * p_reads_buffer;
    P_CordsBuffer &     buffer3 = * p_cords_buffer;
    P_BamLinkBuffer &   buffer4 = * p_bam_link_buffer;
    int n_in = std::min(buffer2.size(assign_it2, buffer2.inIt()), p_parms.thd_assign_num);
    int n_out = std::min(buffer3.size() - buffer3.usedSize(), p_parms.thd_assign_num);
    int n = std::min(n_in, n_out);
    for (int i = 0; i < n; i++)
    {
        appendValue(tasks[thread_id].p_ins, assign_it2); //assign reads
        appendValue(tasks[thread_id].p_outs, assign_it3); //assign cords and bam_Link
        //buffer2.setProtected(assign_it2);
        //Warn here protected must be set before call nextIn to avoid print new added empty buffer.
        buffer3.setProtected(assign_it3); //protect cords
        buffer4.setProtected(assign_it3); //protect bam_link
        buffer3.nextIn();
        buffer4.nextIn();
        buffer2.nextIt(assign_it2);
        buffer3.nextIt(assign_it3);
        dout << "acrs2" << assign_it2 << buffer2.inIt() << buffer3.isEmpty() << "\n";
    }
    int return_val = length(tasks[thread_id].p_ins) ? 0 : 1;
    dout << "acrs1" << length(tasks[thread_id].p_ins) << n << n_in << n_out << "it" << buffer2.inIt() << assign_it2 << buffer2.inIt() << buffer2.outIt() << "\n";
    return return_val;
}

/*----------  Class P_Parms  ----------*/

P_Parms::P_Parms():
    thd_fetch_num(1),
    thd_buffer_block_size(5),
    thd_assign_num(1),
    thd_print_num(1)
    {}
void P_Parms::printParms()
{
    dout << "thd_fetch_num" << thd_fetch_num << "\n";
    dout << "thd_buffer_block_size" << thd_buffer_block_size << "\n";
    dout << "thd_assign_num" << thd_assign_num << "\n";
    dout << "thd_print_num" << thd_print_num << "\n";
}
/*----------  Threads functions  ----------*/
/*
 * Free critical variables 
 */
int freeCurrentTask_(P_Tasks & p_tasks, P_Parms & p_parms, int thread_id, int f_error)
{
    if (f_error)
    {
        p_tasks.setTaskAllEnd();
        p_tasks.setTaskEnd(thread_id);
        //dout << "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx1" << "\n";
    }
    if (p_tasks.isTaskFetchReads(thread_id))
    {
        p_tasks.setTaskWait(thread_id);
        atomicCas(p_tasks.sgn_fetch, 1, 0);
        //dout << "tasks=============================" << p_tasks.tasks[thread_id].task_type << "\n";
    }
    else if (p_tasks.isTaskCalRecords(thread_id))
    {
        clear(p_tasks.tasks[thread_id].p_ins);
        clear(p_tasks.tasks[thread_id].p_outs);
        p_tasks.setTaskWait(thread_id);
        //dout << "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx2" << "\n";
    }
    else if (p_tasks.isTaskPrintResults(thread_id))
    {
        p_tasks.setTaskWait(thread_id);
        atomicCas(p_tasks.sgn_print, 1, 0);
        //dout << "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx3" << "\n";
    }
    else
    {
        //dout << "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx4" << p_tasks.tasks[thread_id].task_type << "\n";
        //none
    }
    return 0;
}
/*
 * Called after calling function @freeCurrentTask_
 */
int requestNewTask_(P_Tasks & p_tasks, P_Parms & p_parms, int thread_id, int f_error)
{
    if (p_tasks.isTaskAllEnd())  
    {
        p_tasks.setTaskEnd(thread_id);
    }
    else if (!atomicCas(p_tasks.sgn_fetch_end, 0, p_tasks.sgn_fetch_end) && 
             !p_tasks.p_reads_buffer->isFull() && 
             !atomicCas(p_tasks.sgn_fetch, 0, 1))
    {
        p_tasks.setTaskFetchReads(thread_id);
    dout << "fetch_t..........." << thread_id << p_tasks.sgn_fetch_end << p_tasks.p_reads_buffer->inIt() << "f" << p_tasks.sgn_fetch << p_tasks.isTaskFetchReads(0) << p_tasks.isTaskFetchReads(1) << "\n";
    }
    else if (!p_tasks.p_cords_buffer-> isEmpty() &&
            !atomicCas(p_tasks.sgn_print, 0, 1))
    {
        p_tasks.setTaskPrintResults(thread_id);
    }
    else if (!p_tasks.assignCalRecords(p_parms, thread_id))
    {
        dout << "cal" << thread_id << "\n";
        p_tasks.setTaskCalRecords(thread_id);
    }
    else
    {
        if (p_tasks.sgn_fetch_end == 1 && 
            p_tasks.p_reads_ids_buffer-> isEmpty() && 
            p_tasks.p_cords_buffer -> isEmpty()) // !!!Add all end conditions here !!!
        {
            dout << "$$$$$$$$$$$$$$$$$$$$$$" << "\n";
            p_tasks.setTaskEnd(thread_id);
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

int p_FetchReads(P_Tasks & p_tasks, P_Parms & p_parms, int thread_id)
{
    P_ReadsIdsBuffer & buffer1 = *p_tasks.p_reads_ids_buffer;
    P_ReadsBuffer & buffer2 = *p_tasks.p_reads_buffer;
    if (p_tasks.paths2_it >= length(p_tasks.paths2))
    {
        //dout << "7&&&&&&&&&&&&&&&&&" << "\n";
        p_tasks.sgn_fetch_end = 1;
        return 0;
    }
    std::string file_name = p_tasks.paths2[p_tasks.paths2_it];
    //dout << "path2" << buffer1.size() << buffer1.inIt() << buffer1.outIt() << p_tasks.paths2_it << length(p_tasks.paths2) << thread_id << file_name << buffer1.isFull() << "\n";
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

        std::ostringstream ss;
        dout << "path2=============" << buffer1.inIt() << buffer1.physicalSize() << thread_id << p_tasks.tasks[0].task_type << p_tasks.tasks[1].task_type << "\n";
        //std::cout << ss.str();

        readRecords(buffer1.in(), buffer2.in(), p_tasks.fin, p_parms.thd_buffer_block_size);
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
        dout << "path3---------------" << buffer1.outIt() << buffer1.inIt() << thread_id << buffer1.isFull() << buffer1.usedSize() << "\n";
        if (atEnd(p_tasks.fin))
        {
            close(p_tasks.fin);
            p_tasks.paths2_it++;
            p_tasks.f_fin_open = 0;
        }
    }
    return 0;
}
/*
 * Calculate records
 */
int p_CalRecords(P_Tasks & p_tasks, P_Parms & p_parms, P_Mapper & p_mapper int thread_id)
{
    P_ReadsIdsBuffer & buffer1 = * p_tasks.p_reads_ids_buffer;
    P_ReadsBuffer & buffer2 = * p_tasks.p_reads_buffer;
    P_CordsBuffer & buffer3 = * p_tasks.p_cords_buffer;
    P_BamLinkBuffer & buffer4 = * p_tasks.p_bam_link_buffer;
    if (buffer1.isEmpty() || buffer3.isFull())
    {
        return 0;   
    }
    //p_tasks.p_reads_buffer -> nextOut();
    //p_tasks.p_reads_ids_buffer -> nextOut();
    for (int i = 0; i < p_tasks.getTaskInBufferLen(thread_id); i++)
    {
        //add cal function here
        int i_out = p_tasks.getTaskOutBufferPtr(thread_id, i);
        int i_in = p_tasks.getTaskInBufferPtr(thread_id, i);
        p_mapper.calRecords();
        dout << "p_calreocords" << i << length(p_tasks.tasks[thread_id].p_ins) << buffer1[p_tasks.tasks[thread_id].p_ins[i]][0] << "\n";
        buffer3.unsetProtected(i_out);
        buffer4.unsetProtected(i_out);
    }

    return 0;
}

int p_PrintResults(P_Tasks & p_tasks, P_Parms p_parms, P_Mapper & p_mapper, int thread_id)
{   
    P_ReadsIdsBuffer & buffer1 = * p_tasks.p_reads_ids_buffer;
    P_ReadsBuffer & buffer2 = * p_tasks.p_reads_buffer;
    P_CordsBuffer & buffer3 = * p_tasks.p_cords_buffer;
    P_BamLinkBuffer & buffer4 = * p_tasks.p_bam_link_buffer;
    if (buffer3.isEmpty())
    {
        return 0;
    }
    dout << "print_tie" << buffer3.outIt() << !buffer3.isProtected(buffer3.outIt()) << !buffer3.isEmpty() <<  "\n";
    for (int i = 0; i < p_parms.thd_print_num && !buffer3.isProtected(buffer3.outIt()) 
        && !buffer3.isEmpty(); i++)
    {
        //add print_function here
        //dout << "p_PrintResults <<<<<<" << buffer1[0][0] << buffer1.outIt() << buffer1.inIt() << "\n";
        p_mapper.PrintResults();
        //reads(ids) and cords(bam) buffer are not released until the following sentences.
        buffer1.nextOut();
        buffer2.nextOut();
        buffer3.nextOut();
        buffer4.nextOut();
    }
    return 0; 
}

int p_ThreadProcess(P_Tasks & p_tasks, P_Parms p_parms, P_Mapper & p_mapper, int thread_id)
{
        //std::cout << "p1" << "\n";
    int f_error = 0;
    while (true)
    {
        if (p_RequestTask(p_tasks, p_parms, thread_id, f_error))
        {   
            continue; //continue request if failed in requesting task
        }
        if (p_tasks.isTaskEnd(thread_id))
        {
            break;
        }
        else if(p_tasks.isTaskFetchReads(thread_id)) 
        {
            f_error = p_FetchReads(p_tasks, p_parms, thread_id);
            dout << "fe1" << f_error << "\n";
        }
        else if (p_tasks.isTaskCalRecords(thread_id))
        {
            f_error = p_CalRecords(p_tasks, p_parms, thread_id);
            dout << "fe2" << f_error << "\n";
        }
        
        else if (p_tasks.isTaskPrintResults(thread_id))
        {
            f_error = p_PrintResults(p_tasks, p_parms, thread_id);
            dout << "fe3" << f_error << "\n";
        }
    } 
    return 0;
}
