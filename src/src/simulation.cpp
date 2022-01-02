/* Portions Copyright 2019-2021 Xuesong Zhou and Peiheng Li, Cafer Avci
 
 * If you help write or modify the code, please also list your names here.
 * The reason of having Copyright info here is to ensure all the modified version, as a whole, under the GPL
 * and further prevent a violation of the GPL.
 *
 * More about "How to use GNU licenses for your own software"
 * http://www.gnu.org/licenses/gpl-howto.html
 */

// Peiheng, 02/03/21, remove them later after adopting better casting
#pragma warning(disable : 4305 4267 4018)
// stop warning: "conversion from 'int' to 'float', possible loss of data"
#pragma warning(disable: 4244)

#ifdef _WIN32
#include "pch.h"
#endif

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <string>
#include <cstring>
#include <cstdio>
#include <ctime>
#include <cmath>
#include <algorithm>
#include <functional>
#include <stack>
#include <list>
#include <vector>
#include <map>
#include <omp.h>
#include "config.h"
#include "utils.h"


using std::max;
using std::min;
using std::cout;
using std::endl;
using std::string;
using std::vector;
using std::map;
using std::ifstream;
using std::ofstream;
using std::istringstream;

#include "DTA.h"

void Assignment::AllocateLinkMemory4Simulation()
{
    g_number_of_simulation_intervals = (g_LoadingEndTimeInMin - g_LoadingStartTimeInMin + _simulation_discharge_period_in_min) * 60 / number_of_seconds_per_interval + 2;

    g_number_of_intervals_in_sec = (g_LoadingEndTimeInMin - g_LoadingStartTimeInMin + _simulation_discharge_period_in_min) * 60;

    dtalog.output() << "LoadingStartTimeInMin = " << g_LoadingStartTimeInMin << endl;
    dtalog.output() << "g_LoadingStartTimeInMin = " << g_LoadingEndTimeInMin << endl;
    dtalog.output() << "number_of_simulation_intervals = " << g_number_of_simulation_intervals << endl;
    dtalog.output() << "number_of_simu intervals in sec = " << g_number_of_intervals_in_sec << endl;

    g_number_of_loading_intervals_in_sec = (g_LoadingEndTimeInMin - g_LoadingStartTimeInMin) * 60;

    g_number_of_intervals_in_min = (int)(g_number_of_simulation_intervals / number_of_simu_intervals_in_min + 1);
    // add + 120 as a buffer
    g_number_of_in_memory_simulation_intervals = g_number_of_simulation_intervals;

    dtalog.output() << "allocate 2D dynamic memory LinkOutFlowCapacity..." << endl;

    m_LinkOutFlowCapacity = Allocate2DDynamicArray <float>(g_number_of_links, g_number_of_intervals_in_sec);  //1
    m_LinkOutFlowState = Allocate2DDynamicArray <int>(g_number_of_links, g_number_of_intervals_in_sec);  //1
    // discharge rate per simulation time interval
    dtalog.output() << "allocate 2D dynamic memory m_LinkCumulativeArrivalVector..." << endl;
    m_LinkCumulativeArrivalVector = Allocate2DDynamicArray <float>(g_number_of_links, g_number_of_intervals_in_min);  //2

    dtalog.output() << "allocate 2D dynamic memory m_LinkCumulativeDepartureVector..." << endl;
    m_LinkCumulativeDepartureVector = Allocate2DDynamicArray <float>(g_number_of_links, g_number_of_intervals_in_min);  //3

    m_LinkCACount = Allocate1DDynamicArray <float>(g_number_of_links);
    m_LinkCDCount = Allocate1DDynamicArray <float>(g_number_of_links);

    dtalog.output() << "allocate 2D dynamic memory m_LinkTDWaitingTime..." << endl;
    m_LinkTDWaitingTime = Allocate2DDynamicArray <float>(g_number_of_links, g_number_of_intervals_in_min); //5

    m_LinkTotalWaitingTimeVector.clear();
    for (int i = 0; i < g_number_of_links; ++i)
    {
        m_LinkTotalWaitingTimeVector.push_back(0.0);
    }

    dtalog.output() << "initializing time dependent capacity data..." << endl;


#pragma omp parallel for
    for (int i = 0; i < g_number_of_links; ++i)
    {
        float cap_count = 0;
        float discharge_rate_per_sec = g_link_vector[i].lane_capacity * g_link_vector[i].number_of_lanes / 3600.0;
        float discharge_rate_after_loading = 10 * discharge_rate_per_sec;  ///to collect travel time statistics * 10 times of capacity to discharge all flow */ at the end of simulation time interval

        if (i == 24)
        {
            int debug = 1;
        }

        unsigned int RandomSeed = 101;
        float residual;
        float random_ratio = 0;

        for (int t = 0; t < g_number_of_intervals_in_sec; ++t)
        {
            float OutFlowRate = 0;
            if (t >= g_number_of_loading_intervals_in_sec)
                OutFlowRate = discharge_rate_after_loading;
            /* 10 times of capacity to discharge all flow */
            else
                OutFlowRate = discharge_rate_per_sec;

            // free flow state
            //cell based simulation 

               m_LinkOutFlowCapacity[i][t] = discharge_rate_after_loading;

                if (g_link_vector[i].cell_type >= 0)  // DTA micro cell based simulation
                    m_LinkOutFlowCapacity[i][t] = 1 * g_link_vector[i].number_of_lanes;
                else
                {
                    residual = OutFlowRate - (int)(OutFlowRate);
                    //RandomSeed is automatically updated.
                    RandomSeed = (LCG_a * RandomSeed + LCG_c) % LCG_M;
                    random_ratio = float(RandomSeed) / LCG_M;

                    if (random_ratio < residual)
                        m_LinkOutFlowCapacity[i][t] = (int)(OutFlowRate)+1;
                    else
                        m_LinkOutFlowCapacity[i][t] = (int)(OutFlowRate);

                    if (t < g_number_of_loading_intervals_in_sec)
                        cap_count += m_LinkOutFlowCapacity[i][t];

                }
            // 1800 rule
            //if(t%2==0)
            //    m_LinkOutFlowCapacity[i][t] = 1*g_link_vector[i].number_of_lanes;
            //if (t % 2 == 1)
            //    m_LinkOutFlowCapacity[i][t] = 0;



        }

        for (int t = 0; t < g_number_of_intervals_in_min; ++t)
        {        // convert per hour capacity to per second and then to per simulation interval
            m_LinkCumulativeArrivalVector[i][t] = 0;
            m_LinkCumulativeDepartureVector[i][t] = 0;
            m_LinkTDWaitingTime[i][t] = 0;
            m_LinkOutFlowState[i][t] = 1;
        }


    }

    // for each link with  for link type code is 's'
    //reset the time-dependent capacity to zero
    //go through the records defined in timing_arc file
    //only enable m_LinkOutFlowCapacity[l][t] for the timestamp in the time window per time interval and reset the link TD travel time.


    for (unsigned l = 0; l < g_link_vector.size(); ++l)
    {
        if (l == 40949)
        {
            int i_bebug = 1;
        }
        if (g_link_vector[l].timing_arc_flag == true)  // with timing data
        {

            if (g_link_vector[l].timing_arc_flag)
            {
                // reset for signalized links (not freeway links as type code != 'f' for the case of freeway workzones)
                // only for the loading period
                for (int t = 0; t < g_number_of_loading_intervals_in_sec; ++t)
                {
                    m_LinkOutFlowCapacity[l][t] = 0;
                    m_LinkOutFlowState[l][t] = 0;

                }
            }

            int number_of_cycles = (g_LoadingEndTimeInMin - g_LoadingStartTimeInMin) * 60 / max(1, g_link_vector[l].VDF_period[0].cycle_link_distance_in_km);  // unit: seconds;

            int cycle_link_distance_in_km = g_link_vector[l].VDF_period[0].cycle_link_distance_in_km;
            int start_green_time = g_link_vector[l].VDF_period[0].start_green_time;
            int end_green_time = g_link_vector[l].VDF_period[0].end_green_time;

            if (end_green_time < start_green_time)
            {
                end_green_time += cycle_link_distance_in_km;  // consider a looped end green time notation, e.g. 60-10 for cl = 100, then end green time should be 110. 
            }
            dtalog.output() << "signal timing data: link: cycle_link_distance_in_km = " << cycle_link_distance_in_km <<
                ",start_green_time = " << start_green_time <<
                ",end_green_time = " << end_green_time <<
                endl;

            for (int cycle_no = 0; cycle_no < number_of_cycles; ++cycle_no)
            {
                int count = 0;

                // relative time horizon
                for (int t = cycle_no * cycle_link_distance_in_km + start_green_time; t <= cycle_no * cycle_link_distance_in_km + end_green_time; t += 1)
                {
                    // activate capacity for this time duration
                    m_LinkOutFlowCapacity[l][t] = g_link_vector[l].saturation_flow_rate * g_link_vector[l].number_of_lanes / 3600.0;

                    if (t % 2 == 0)
                        m_LinkOutFlowCapacity[l][t] = 1 * g_link_vector[l].number_of_lanes;
                    if (t % 2 == 1)
                        m_LinkOutFlowCapacity[l][t] = 0;

                    m_LinkOutFlowState[l][t] = 1;

                    //residual = m_LinkOutFlowCapacity[l][t] - (int)(m_LinkOutFlowCapacity[l][t]);
                    ////RandomSeed is automatically updated.
                    //RandomSeed = (LCG_a * RandomSeed + LCG_c) % LCG_M;
                    //random_ratio = float(RandomSeed) / LCG_M;

                    //if (random_ratio < residual)
                    //    m_LinkOutFlowCapacity[l][t] = (int)(m_LinkOutFlowCapacity[l][t]) + 1;
                    //else
                    //    m_LinkOutFlowCapacity[l][t] = (int)(m_LinkOutFlowCapacity[l][t]);
                }
            }
        } if (g_link_vector[l].m_link_pedefined_capacity_map.size() > 0)
        {


            for (int t = 0; t < g_number_of_loading_intervals_in_sec; ++t)
            {
                if(g_link_vector[l].m_link_pedefined_capacity_map.find(t) != g_link_vector[l].m_link_pedefined_capacity_map.end())
                { 
                        // predefined outlflow capacity 
                    m_LinkOutFlowCapacity[l][t] = g_link_vector[l].m_link_pedefined_capacity_map[t];
                    m_LinkOutFlowState[l][t] = g_link_vector[l].m_link_pedefined_capacity_map[t];
                }
            }
        }
        
        //end f
    }

    dtalog.output() << "End of initializing time dependent capacity data." << endl;
}

void Assignment::DeallocateLinkMemory4Simulation()
{
    // g_fout << "deallocate 2D dynamic memory m_LinkOutFlowCapacity..." << endl;
    //if(m_LinkOutFlowCapacity)
    //    Deallocate2DDynamicArray(m_LinkOutFlowCapacity , g_number_of_links, g_number_of_intervals_in_sec);  //1

   // memory leak here: Xuesong: 11/20/2021                                                                                                         // 
// //if (m_LinkOutFlowState)
    //    Deallocate2DDynamicArray(m_LinkOutFlowState, g_number_of_links, g_number_of_intervals_in_sec);  //1
    // g_fout << "deallocate 2D dynamic memory m_LinkCumulativeArrivalVector..." << endl;
    if (m_LinkCumulativeArrivalVector)
        Deallocate2DDynamicArray(m_LinkCumulativeArrivalVector, g_number_of_links, g_number_of_intervals_in_min);  //2
    // g_fout << "deallocate 2D dynamic memory m_LinkCumulativeDepartureVector..." << endl;
    if (m_LinkCumulativeDepartureVector)
        Deallocate2DDynamicArray(m_LinkCumulativeDepartureVector, g_number_of_links, g_number_of_intervals_in_min); //3

    if (m_LinkCACount)
        Deallocate1DDynamicArray(m_LinkCACount, g_number_of_links);

    if (m_LinkCDCount)
        Deallocate1DDynamicArray(m_LinkCDCount, g_number_of_links);

    if (m_LinkTDWaitingTime)
        Deallocate2DDynamicArray(m_LinkTDWaitingTime, g_number_of_links, g_number_of_intervals_in_min); //4
}

void Assignment::STTrafficSimulation()
{

    ofstream simu_log_file;
    simu_log_file.open("model_simu_log.txt");
    simu_log_file << "simulation log file.\n";

    //given p_agent->path_link_seq_no_vector path link sequence no for each agent
    int TotalCumulative_Arrival_Count = 0;
    int TotalCumulative_Departure_Count = 0;
    double TotalVehicleMileTraveled = 0;
    double TotalVehicleHourTraveled = 0;

    clock_t start_t;
    start_t = clock();

    AllocateLinkMemory4Simulation();

    int agent_type_size = g_AgentTypeVector.size();
    int zone_size = g_zone_vector.size();
    int demand_period_size = g_DemandPeriodVector.size();

    CColumnVector* p_column_pool;
    float path_toll = 0;
    float path_distance = 0;
    float path_travel_time = 0;
    float time_stamp = 0;

    int trace_agent_id = 1105;  //default can be -1
    int trace_link_no = -1;  //default can be -1
    int trace_node_no = -100;  //default can be -1

    std::map<int, CColumnPath>::iterator it, it_begin, it_end;

    for (int orig = 0; orig < zone_size; ++orig)
    {
        if (orig % 100 == 0)
            dtalog.output() << "generating " << g_agent_simu_vector.size() / 1000 << " K agents for " << orig << "  zones " << endl;

        for (int at = 0; at < agent_type_size; ++at)
        {
            for (int dest = 0; dest < zone_size; ++dest)
            {
                for (int tau = 0; tau < demand_period_size; ++tau)
                {
                    p_column_pool = &(assignment.g_column_pool[orig][dest][at][tau]);
                    if (p_column_pool->od_volume > 0)
                    {
                        // scan through the map with different node sum for different continuous paths
                        it_begin = p_column_pool->path_node_sequence_map.begin();
                        it_end = p_column_pool->path_node_sequence_map.end();

                        for (it = it_begin; it != it_end; ++it)
                        {
                            path_toll = 0;
                            path_distance = 0;
                            path_travel_time = 0;

                            int VehicleSize = (it->second.path_volume + 0.5);

                            for (int v = 0; v < VehicleSize; ++v)
                            {

                                int slot_no = assignment.g_DemandPeriodVector[tau].get_time_slot_no();

                                if (it->second.path_volume < 1)
                                    time_stamp = (slot_no - assignment.g_DemandPeriodVector[tau].starting_time_slot_no) * MIN_PER_TIMESLOT + g_agent_simu_vector.size() % 15 / 15.0 * 1.0 * MIN_PER_TIMESLOT;
                                else
                                    time_stamp = (slot_no - assignment.g_DemandPeriodVector[tau].starting_time_slot_no) * MIN_PER_TIMESLOT + v * 1.0 / it->second.path_volume * MIN_PER_TIMESLOT;

                                if (it->second.m_link_size == 0)   // only load agents with physical path
                                    continue;

                                CAgent_Simu* pAgent = new CAgent_Simu();

                                if (pAgent->agent_id == 45)
                                {
                                    int idebug = 1;
                                }
                                // for future use of column pool
                                pAgent->at = at;
                                pAgent->dest = dest;
                                pAgent->tau = tau;

                                pAgent->PCE_unit_size = max(1, (int)(assignment.g_AgentTypeVector[at].PCE + 0.5));  // convert a possible floating point pce to an integer value for simulation
                                pAgent->time_headway = (int)(assignment.g_AgentTypeVector[at].time_headway_in_sec / number_of_seconds_per_interval + 0.5);
                                pAgent->agent_id = g_agent_simu_vector.size();
                                pAgent->departure_time_in_min = time_stamp;

                                pAgent->path_travel_time_in_min = g_LoadingEndTimeInMin - pAgent->departure_time_in_min;  // by default
                                it->second.agent_simu_id_vector.push_back(pAgent->agent_id);

                                int simulation_time_intervalNo = (int)(pAgent->departure_time_in_min * 60 / number_of_seconds_per_interval);

                                if (simulation_time_intervalNo < 0)  // reset to a zero index departure time interval
                                {
                                    simulation_time_intervalNo = 0;
                                    int idebug = 1;
                                }

                                if (simulation_time_intervalNo > 1440 * 60 * 4)
                                {
                                    int idebug = 1;
                                }


                                g_AgentTDListMap[simulation_time_intervalNo].m_AgentIDVector.push_back(pAgent->agent_id);

                                for (int nl = 0; nl < it->second.m_link_size; ++nl)  // arc a
                                {
                                    int link_seq_no = it->second.path_link_vector[nl];
                                    pAgent->path_link_seq_no_vector.push_back(link_seq_no);
                                }
                                pAgent->AllocateMemory();

                                int FirstLink = pAgent->path_link_seq_no_vector[0];

                                pAgent->m_Veh_LinkArrivalTime_in_simu_interval[0] = simulation_time_intervalNo;
                                pAgent->m_Veh_LinkDepartureTime_in_simu_interval[0] = pAgent->m_Veh_LinkArrivalTime_in_simu_interval[0] + (int)(g_link_vector[FirstLink].free_flow_travel_time_in_min * number_of_simu_intervals_in_min);

                                g_agent_simu_vector.push_back(pAgent);
                            }
                        }
                    }
                }
            }
        }
    }

    dtalog.output() << "number of simulation zones:" << zone_size << endl;

    int current_active_agent_id = 0;
    // the number of threads is redifined.
    int number_of_threads = omp_get_max_threads();

    int number_of_in_min_for_RTSP = 5;

    int link_size = g_link_vector.size();

    // initialize CA and CD count for each link

    bool cell_based_simulation_mode = false;
    for (int li = 0; li < link_size; ++li)
    {
        CLink* pLink = &(g_link_vector[li]);

        if (pLink->cell_type >= 0)  // activate cell based simulation mode
            cell_based_simulation_mode = true;

        m_LinkCACount[li] = 0;
        m_LinkCDCount[li] = 0;

        g_link_vector[li].win_count = 0;
        g_link_vector[li].lose_count = 0;
        g_link_vector[li].time_to_be_released = -1;

    }


    bool bRealTimeInformationActivated = false;
    // first loop for time t
    for (int t = 0; t < g_number_of_simulation_intervals; ++t)
    {

        if (t % number_of_simu_intervals_in_min == 0)  // every min
        {
            int time_in_min = t / number_of_simu_intervals_in_min;
            for (int li = 0; li < link_size; ++li)
            {
                CLink* pLink = &(g_link_vector[li]);
                {
                    m_LinkCumulativeArrivalVector[li][time_in_min] = m_LinkCACount[li];
                    m_LinkCumulativeDepartureVector[li][time_in_min] = m_LinkCDCount[li];
                }
            }
        }

        if (t % (int(number_of_in_min_for_RTSP * 60 / number_of_seconds_per_interval)) == 0)
        {
            int slot_no = (t * number_of_seconds_per_interval / 60 + g_LoadingStartTimeInMin) / 15;
            if (RTSP_RealTimeShortestPathFinding(slot_no, t) == true)
                bRealTimeInformationActivated = true;
        }

        int number_of_simu_interval_per_min = 60 / number_of_seconds_per_interval;
        double network_wide_speed = 60;
        double network_wide_travel_time = 0;

        if (TotalVehicleHourTraveled > 0.001)
        {
            network_wide_speed = TotalVehicleMileTraveled / max(0.0001, TotalVehicleHourTraveled);
            network_wide_travel_time = TotalVehicleHourTraveled / max(1, TotalCumulative_Departure_Count) * 60.0; // from hour to min
        }


        if (t % (number_of_simu_interval_per_min * 10) == 0)
        {
            dtalog.output() << "simu time= " << t / number_of_simu_interval_per_min << " min, CA = " << TotalCumulative_Arrival_Count << " CD=" << TotalCumulative_Departure_Count
                << ", speed=" << network_wide_speed << ", travel time = " << network_wide_travel_time
                << endl;

            simu_log_file << "simu time= " << t / number_of_simu_interval_per_min << " min, CA = " << TotalCumulative_Arrival_Count << " CD=" << TotalCumulative_Departure_Count
                << ", speed=" << network_wide_speed << ", travel time = " << network_wide_travel_time
                << endl;
        }

        if (g_AgentTDListMap.find(t) != g_AgentTDListMap.end())
        {
            for (int Agent_v = 0; Agent_v < g_AgentTDListMap[t].m_AgentIDVector.size(); ++Agent_v)
            {
                int agent_id = g_AgentTDListMap[t].m_AgentIDVector[Agent_v];

                CAgent_Simu* p_agent = g_agent_simu_vector[agent_id];
                p_agent->m_bGenereated = true;
                int FirstLink = p_agent->path_link_seq_no_vector[0];
                m_LinkCACount[FirstLink] += 1;

                g_link_vector[FirstLink].EntranceQueue.push_back(p_agent->agent_id);


                if (trace_agent_id == agent_id)
                {
                    simu_log_file << "trace tag 1: simu time= " << t / number_of_simu_interval_per_min << " min, , traced vehicle enters the network = "  << endl;
                }


                TotalCumulative_Arrival_Count++;
            }
        }

#pragma omp parallel for    // parallel computing for each link
        for (int li = 0; li < link_size; ++li)
        {
            CLink* pLink = &(g_link_vector[li]);

            // if there are Agents in the entrance queue
            while (pLink->EntranceQueue.size() > 0)
            {
                int agent_id = pLink->EntranceQueue.front();

                pLink->EntranceQueue.pop_front();
                pLink->ExitQueue.push_back(agent_id);
                CAgent_Simu* p_agent = g_agent_simu_vector[agent_id];
                if (trace_agent_id == agent_id)
                {
                    simu_log_file << "trace tag 2: simu time interval = " << t << " min, , traced vehicle moves from entrance queue to exit queue on link = "
                        << g_node_vector[pLink->from_node_seq_no].node_id << " -> " << g_node_vector[pLink->to_node_seq_no].node_id <<
                        " on its link seq.no " << p_agent->m_current_link_seq_no << endl;

                    trace_link_no = li;
                    trace_node_no = g_node_vector[pLink->to_node_seq_no].node_seq_no;
                }

                int arrival_time_in_simu_interval = p_agent->m_Veh_LinkArrivalTime_in_simu_interval[p_agent->m_current_link_seq_no];
                int link_travel_time_in_simu_interavls = max(1, g_link_vector[li].free_flow_travel_time_in_min * number_of_simu_intervals_in_min);
                p_agent->m_Veh_LinkDepartureTime_in_simu_interval[p_agent->m_current_link_seq_no] = arrival_time_in_simu_interval + g_link_vector[li].free_flow_travel_time_in_min * number_of_simu_intervals_in_min;

                //dtalog.output() << "reserve TD at time t = " << t << " on link" << g_node_vector[pLink->from_node_seq_no].node_id << " -> " << g_node_vector[pLink->to_node_seq_no].node_id <<
                //    " for agent = " << agent_id << " on link seq. no " << p_agent->m_current_link_seq_no << " with travel time in simu intervas = " << link_travel_time_in_simu_interavls << endl;

            }
        }

        /// for each link, for all vehicles in the queue 
        ///                         // back trace the resourece use vector right after the vehicle moves to the next link
        for (int li = 0; li < link_size; ++li)
        {
            CLink* pLink = &(g_link_vector[li]);
            int queue_position_no = 0;
            for (auto it = pLink->ExitQueue.begin(); it != pLink->ExitQueue.end(); ++it)
            {
                int agent_id = (*it);
                CAgent_Simu* p_agent = g_agent_simu_vector[agent_id];

                if (trace_agent_id == agent_id)
                {
                    simu_log_file << "trace tag 2.5: simu time= " << t / number_of_simu_interval_per_min << " min, , traced vehicle is waiting at the exit queue on link = "
                        << g_node_vector[pLink->from_node_seq_no].node_id << " -> " << g_node_vector[pLink->to_node_seq_no].node_id <<
                        " on its link seq.no " << p_agent->m_current_link_seq_no <<
                        " queue pos.no " << queue_position_no <<
                        " of queue size = " << pLink->ExitQueue.size() << endl;

                    trace_link_no = li;
                }

                queue_position_no++;

                if (p_agent->PCE_unit_size >= 2)
                {
                    for (int l_backtrace = 1; l_backtrace < p_agent->PCE_unit_size; l_backtrace++)
                    {
                        int local_l_index = p_agent->m_current_link_seq_no - l_backtrace;
                        if (local_l_index >= 0)
                        {
                            int current_link_seq_no = p_agent->path_link_seq_no_vector[local_l_index];
                            CLink* pCurrentLink = &(g_link_vector[current_link_seq_no]);
                            pCurrentLink->time_to_be_released = t + p_agent->time_headway;
                            // big truck/bus, 
                            //backtrace the previous l - k links, k = 1, 2, K, and set release time
                              //  K as the number of units of truck
                        }

                    }
                }
            }
        } // end of link 

                /// 

        /// 
        /// 
        /// </summary>
        int node_size = g_node_vector.size();

#pragma omp parallel for  // parallel computing for each node
        for (int node = 0; node < node_size; ++node)  // for each node
        {
            if (node== trace_node_no)
            {
                simu_log_file << "trace tag 3F_1: simu time int = " << t  << ", , check node incoming mode "
                    << g_node_vector[node].node_id <<  endl;
            }

            bool node_resource_competing_mode = false;
            int incoming_link_request_count = 0;
            int incoming_link_index_FIFO = -1;
            int debug_node_resource_competing_mode = 0;
            if (cell_based_simulation_mode)
            {
                g_node_vector[node].next_link_for_resource_request.clear();


                CLink* pLink;
                CLink* pNextLink;
                for (int i = 0; i < g_node_vector[node].m_incoming_link_seq_no_vector.size(); ++i)
                {
                    int link = g_node_vector[node].m_incoming_link_seq_no_vector[i];  // we will start with different first link from the incoming link list,
                    // equal change, regardless of # of lanes or main line vs. ramp, but one can use service arc, to control the effective capacity rates, e.g. through a metered ramp, to
                    // allow mainline to use the remaining flow
                    pLink = &(g_link_vector[link]);

                    if (node == trace_node_no)
                    {
                        simu_log_file << "trace tag 3F_2: simu time int = " << t << ", , check node incoming mode "
                            << g_node_vector[node].node_id << ",i= " << i << ", link id"
                            << g_node_vector[pLink->from_node_seq_no].node_id << " -> " << g_node_vector[pLink->to_node_seq_no].node_id

                            << endl;
                    }



                    if (pLink->ExitQueue.size() >= 1)
                    {
                        int agent_id = pLink->ExitQueue.front();
                        CAgent_Simu* p_agent = g_agent_simu_vector[agent_id];
                        int next_link_seq_no = p_agent->path_link_seq_no_vector[p_agent->m_current_link_seq_no + 1];
                         pNextLink = &(g_link_vector[next_link_seq_no]);

                         if (trace_agent_id == agent_id)
                         {
                             simu_log_file << "trace tag 3: simu time= " << t / number_of_simu_interval_per_min << " min, , traced vehicle wants to transfer from link = "
                                 << g_node_vector[pLink->from_node_seq_no].node_id << " -> " << g_node_vector[pLink->to_node_seq_no].node_id <<
                                 " on its link seq.no " << p_agent->m_current_link_seq_no << " to the next link. " << endl;
                         }

                        int current_vehicle_count = m_LinkCACount[next_link_seq_no] - m_LinkCDCount[next_link_seq_no];

                        if (trace_agent_id == agent_id)
                        {
                            simu_log_file << "trace tag 3.5: simu time= " << t*1.0 / number_of_simu_interval_per_min << " min, , traced vehicle checks next link = "
                                << g_node_vector[pNextLink->from_node_seq_no].node_id << " -> " << g_node_vector[pNextLink->to_node_seq_no].node_id <<
                                ", cell type = " << pNextLink->cell_type << ", current_vehicle_count = " << current_vehicle_count <<
                                ",pNextLink->spatial_capacity_in_vehicles= " << pNextLink->spatial_capacity_in_vehicles << endl;

                        }

                        if (pNextLink->cell_type >= 0 && current_vehicle_count < pNextLink->spatial_capacity_in_vehicles)  // only apply for cell mode
                        {
                            g_node_vector[node].next_link_for_resource_request[next_link_seq_no] = 1;
                            incoming_link_request_count++;
                        }
                    }
                }


                if (incoming_link_request_count >= 2)
                {
                    if (g_node_vector[node].next_link_for_resource_request.size() >= 1)
                    {
                        //                if(g_node_vector[node].node_id == 2347 && t / 240 >= 18)
                        node_resource_competing_mode = true;

                    }
                }

                if (node == trace_node_no)
                {
                    simu_log_file << "trace tag 3F_2: simu time= " << t * 1.0 / number_of_simu_interval_per_min << " min, , check node incoming mode "
                        << g_node_vector[node].node_id
                        << "next_link_for_resource_request size = " << g_node_vector[node].next_link_for_resource_request.size()
                        << "node_resource_competing_mode = " << node_resource_competing_mode
                        << endl;
                }

                // determine FIFO link with earliest departure time request
                if (node_resource_competing_mode)
                {
                    int earlest_departure_time_interval = 99999999;

                    for (int i = 0; i < g_node_vector[node].m_incoming_link_seq_no_vector.size(); ++i)
                    {
                        int link = g_node_vector[node].m_incoming_link_seq_no_vector[i];  // we will start with different first link from the incoming link list,
                        // equal change, regardless of # of lanes or main line vs. ramp, but one can use service arc, to control the effective capacity rates, e.g. through a metered ramp, to
                        // allow mainline to use the remaining flow
                        pLink = &(g_link_vector[link]);

                        if (trace_link_no == link)
                        {
                                simu_log_file << "trace tag 3.6E: simu time= " << t * 1.0 / number_of_simu_interval_per_min << " min, , traced vehicle uses FIFO rule on link "
                                    << g_node_vector[pLink->from_node_seq_no].node_id << " -> " << g_node_vector[pLink->to_node_seq_no].node_id <<
                                    " queue size = " << pLink->ExitQueue.size() << endl;
                        }
                        if (pLink->cell_type >= 0 && pLink->ExitQueue.size() >= 1)
                        {
                            int agent_id = pLink->ExitQueue.front();
                            CAgent_Simu* p_agent = g_agent_simu_vector[agent_id];
                            if (p_agent->m_Veh_LinkDepartureTime_in_simu_interval[p_agent->m_current_link_seq_no] < earlest_departure_time_interval)
                            {

                                if (trace_agent_id == agent_id)
                                {
                                    simu_log_file << "trace tag 3.6A: simu time= " << t * 1.0 / number_of_simu_interval_per_min << " min, , traced vehicle uses FIFO rule on link "
                                        << g_node_vector[pLink->from_node_seq_no].node_id << " -> " << g_node_vector[pLink->to_node_seq_no].node_id <<
                                        " expectd departure time = " << p_agent->m_Veh_LinkDepartureTime_in_simu_interval[p_agent->m_current_link_seq_no]<< 
                                        endl;
                                }

                                earlest_departure_time_interval = p_agent->m_Veh_LinkDepartureTime_in_simu_interval[p_agent->m_current_link_seq_no];
                                incoming_link_index_FIFO = i;  // i is the link index in the vector of m_incoming_link_seq_no_vector of node
                            }
                            else
                            {
                                if (trace_agent_id == agent_id)
                                {
                                    simu_log_file << "trace tag 3.6B: simu time= " << t * 1.0 / number_of_simu_interval_per_min << " min, , traced vehicle fails to use FIFO rule on link "
                                        << g_node_vector[pLink->from_node_seq_no].node_id << " -> " << g_node_vector[pLink->to_node_seq_no].node_id <<
                                        " expectd departure time = " << p_agent->m_Veh_LinkDepartureTime_in_simu_interval[p_agent->m_current_link_seq_no] <<
                                        " while " <<
                                        " earlest_departure_time_interval = " << earlest_departure_time_interval <<endl;
                                }

                            }

                        }
                    }
                }


            }
            // for each incoming link

            if (node == trace_node_no)
            {
                simu_log_file << "trace tag 3F: simu time= " << t * 1.0 / number_of_simu_interval_per_min << " min, , check node incoming mode "
                    << g_node_vector[node].node_id
                    << "incoming_link_index_FIFO = " << incoming_link_index_FIFO
                    << endl;
            }

            for (int i = 0; i < g_node_vector[node].m_incoming_link_seq_no_vector.size(); ++i)
            {
                int incoming_link_index = i;

                if (incoming_link_index_FIFO >= 0)
                {
                    incoming_link_index = (i + incoming_link_index_FIFO) % (g_node_vector[node].m_incoming_link_seq_no_vector.size());  // cycle loop
                    if (debug_node_resource_competing_mode)
                        dtalog.output() << "FIFO judgement i = " << i << endl;
                }
                else
                {  // random mode
                    incoming_link_index = (i + t) % (g_node_vector[node].m_incoming_link_seq_no_vector.size());
                }


                int link = g_node_vector[node].m_incoming_link_seq_no_vector[incoming_link_index];  // we will start with different first link from the incoming link list,
                // equal change, regardless of # of lanes or main line vs. ramp, but one can use service arc, to control the effective capacity rates, e.g. through a metered ramp, to
                // allow mainline to use the remaining flow
                CLink* pLink = &(g_link_vector[link]);

                if (node == trace_node_no)
                {
                    simu_log_file << "trace tag 3F: simu time= " << t * 1.0 / number_of_simu_interval_per_min << " min, , check node incoming mode "
                        << g_node_vector[node].node_id
                        << ", check link" << g_node_vector[pLink->from_node_seq_no].node_id << " -> " << g_node_vector[pLink->to_node_seq_no].node_id << endl;

                }


                // check if the current link has sufficient capacity
                // most critical and time-consuming task, check link outflow capacity

                int time_in_sec = t * number_of_seconds_per_interval;


                while (m_LinkOutFlowCapacity[link][time_in_sec] >= 1 && pLink->ExitQueue.size() >= 1)
                {
                    int agent_id = pLink->ExitQueue.front();
                    CAgent_Simu* p_agent = g_agent_simu_vector[agent_id];

                    if (p_agent->m_Veh_LinkDepartureTime_in_simu_interval[p_agent->m_current_link_seq_no] > t)
                    {
                        // the future departure time on this link is later than the current time
                        break;
                    }

                    if (p_agent->m_current_link_seq_no == p_agent->path_link_seq_no_vector.size() - 2)  //-2 do not consider virtual destination arc
                    {
                        // end of path
                        pLink->ExitQueue.pop_front();

                        p_agent->m_bCompleteTrip = true;
                        m_LinkCDCount[link] += 1;
                        p_agent->arrival_time_in_min = t * 1.0 / number_of_simu_interval_per_min;
                        p_agent->path_travel_time_in_min = p_agent->arrival_time_in_min - p_agent->departure_time_in_min;

                        if (p_agent->agent_id == 0)
                        {
                            int debug = 1;
                        }
                        p_agent->path_distance = 0;

                        for (int link_s = 0; link_s < p_agent->path_link_seq_no_vector.size(); link_s++)
                        {
                            int link_seq_no = p_agent->path_link_seq_no_vector[link_s];
                            //path_toll += g_link_vector[link_seq_no].VDF_period[p_agent->tau].toll[at];
                            p_agent->path_distance += g_link_vector[link_seq_no].link_distance_in_km;
                        }

#pragma omp critical
                        {
                            TotalCumulative_Departure_Count += 1;
                            TotalVehicleHourTraveled += p_agent->path_travel_time_in_min / 60.0;
                            TotalVehicleMileTraveled += p_agent->path_distance;
                        }

                    }
                    else
                    {
                        //RT re-routing 
                        // test the one to all trees
                        // replease link sequence 
                        // contiue
                        // 
                        // 
                        // not complete the trip. move to the next link's entrance queue
                        int next_link_seq_no = p_agent->path_link_seq_no_vector[p_agent->m_current_link_seq_no + 1];
                        CLink* pNextLink = &(g_link_vector[next_link_seq_no]);

                        if (p_agent->m_current_link_seq_no >= 4)
                        {
                            int debug = 1;
                        }
                        //// spatial queue
                        ////if(pNextLink->link_distance_in_km <=0.008)  // cell based link
                        ////{ 
                        //    pNextLink->spatial_capacity_in_vehicles = 1; // for cell based model
                        ////}

                        if (1/*pNextLink->cell_type >= 0 || pNextLink->traffic_flow_code == 2*/)  // DTA micro cell based simulation
                        {
                            int current_vehicle_count = m_LinkCACount[next_link_seq_no] - m_LinkCDCount[next_link_seq_no];
                            if (current_vehicle_count >= pNextLink->spatial_capacity_in_vehicles || (t < pNextLink->time_to_be_released))
                                // this is an OR condition, related to reaction time tau
                            {
                                // spatical queue in the next link is blocked, break the while loop from here, as a first in first out queue.
                                // g_fout << "spatical queue in the next link is blocked on link seq  " <<  g_node_vector[pNextLink->from_node_seq_no].node_id  << " -> " << g_node_vector[pNextLink->to_node_seq_no].node_id  <<endl;
                                //dtalog.output() << "spatical queue in the next link is blocked at time  = " << t << ",link" << g_node_vector[pLink->from_node_seq_no].node_id << " -> " << g_node_vector[pLink->to_node_seq_no].node_id <<
                                //    ",current_vehicle_count = " << current_vehicle_count << endl;
                                if (node_resource_competing_mode)
                                {
                                    pLink->lose_count++;
                                    if (debug_node_resource_competing_mode)
                                    {
                                        simu_log_file << "lose: the next link is blocked at time  = " << t * 1.0 << ", from link" << g_node_vector[pLink->from_node_seq_no].node_id << " -> " << g_node_vector[pLink->to_node_seq_no].node_id
                                            << " to link" << g_node_vector[pNextLink->from_node_seq_no].node_id << " -> " << g_node_vector[pNextLink->to_node_seq_no].node_id
                                            << " time request " << p_agent->m_Veh_LinkDepartureTime_in_simu_interval[p_agent->m_current_link_seq_no]
                                            << endl;
                                    }
                                }
                                break;  // being blocked
                            }
                            else // free to go if none of the above conditions cannot be met 
                            {
                                if (node_resource_competing_mode)
                                {
                                    pLink->win_count++;
                                    if (debug_node_resource_competing_mode)
                                    {
                                        dtalog.output() << "win:spatical queue in the next link is blocked at time  = " << t << ",link" << g_node_vector[pLink->from_node_seq_no].node_id << " -> " << g_node_vector[pLink->to_node_seq_no].node_id
                                            << " to link" << g_node_vector[pNextLink->from_node_seq_no].node_id << " -> " << g_node_vector[pNextLink->to_node_seq_no].node_id
                                            << " time request " << p_agent->m_Veh_LinkDepartureTime_in_simu_interval[p_agent->m_current_link_seq_no]
                                            << endl;
                                    }
                                }


                            }

                        }

                        //                // kinematic wave
                        //                if (pNextLink->traffic_flow_code == 3)  // to be discussed later with Cafer
                        //                {
                        //                    int lagged_time_stamp = max(0, t - 1 - pNextLink->BWTT_in_simulation_interval);
                        //                    int current_vehicle_count = m_LinkCACount[next_link_seq_no] - m_LinkCDCount[next_link_seq_no];
                        //                    if (current_vehicle_count > pNextLink->spatial_capacity_in_vehicles)
                        //                    {
                        //                        // spatical queue in the next link is blocked, break the while loop from here, as a first in first out queue.
                                                //// g_fout << "spatical queue in the next link is blocked on link seq  " <<  g_node_vector[pNextLink->from_node_seq_no].node_id  << " -> " << g_node_vector[pNextLink->to_node_seq_no].node_id  <<endl;
                        //                        break;
                        //                    }
                        //                }

                        pLink->ExitQueue.pop_front();

                        //if(pLink->ExitQueue.size() !=0 || pLink->EntranceQueue.size() != 0)
                        //{
                        //    simu_log_file << "Error at time = " << t << " , link = "
                        //        << g_node_vector[pLink->from_node_seq_no].node_id << " -> " << g_node_vector[pLink->to_node_seq_no].node_id <<

                        //        << "ExitQueue.size() = " << pLink->ExitQueue.size()
                        //        << "EntranceQueue.size() = " << pLink->EntranceQueue.size() 
                        //        << "spatial queue capacity =  " << pLink->spatial_capacity_in_vehicles << endl;
                        //}

                        pLink->time_to_be_released = t + p_agent->time_headway;

                        // back trace the resourece use vector right after the vehicle moves to the next link
                        if (p_agent->PCE_unit_size >= 2)
                        {
                            for (int l_backtrace = 1; l_backtrace < p_agent->PCE_unit_size; l_backtrace++)
                            {
                                int local_l_index = p_agent->m_current_link_seq_no - l_backtrace;
                                if (local_l_index >= 0)
                                {
                                    int current_link_seq_no = p_agent->path_link_seq_no_vector[local_l_index];
                                    CLink* pCurrentLink = &(g_link_vector[current_link_seq_no]);
                                    pCurrentLink->time_to_be_released = t + p_agent->time_headway;
                                    // big truck/bus, 
                                    //backtrace the previous l - k links, k = 1, 2, K, and set release time
                                      //  K as the number of units of truck
                                }

                            }

                        }

                        if (trace_agent_id == agent_id)
                        {
                            simu_log_file << "trace tag 4: simu time interval = " << t << " , traced vehicle transfers from link = "
                                << g_node_vector[pLink->from_node_seq_no].node_id << " -> " << g_node_vector[pLink->to_node_seq_no].node_id <<
                                "  to next lnk " <<
                                g_node_vector[pNextLink->from_node_seq_no].node_id << " -> " << g_node_vector[pNextLink->to_node_seq_no].node_id << endl;
                        }

                        pNextLink->EntranceQueue.push_back(agent_id);


                        p_agent->m_Veh_LinkDepartureTime_in_simu_interval[p_agent->m_current_link_seq_no] = t;
                        p_agent->m_Veh_LinkArrivalTime_in_simu_interval[p_agent->m_current_link_seq_no + 1] = t;

                        float travel_time_in_sec = (p_agent->m_Veh_LinkDepartureTime_in_simu_interval[p_agent->m_current_link_seq_no] - p_agent->m_Veh_LinkArrivalTime_in_simu_interval[p_agent->m_current_link_seq_no]) * number_of_seconds_per_interval;
                        //for each waited vehicle
                        float waiting_time_in_min = travel_time_in_sec / 60.0 - pLink->free_flow_travel_time_in_min;

                        m_LinkTDWaitingTime[link][p_agent->m_Veh_LinkArrivalTime_in_simu_interval[p_agent->m_current_link_seq_no] / number_of_simu_intervals_in_min] += waiting_time_in_min;
                        m_LinkTotalWaitingTimeVector[link] += waiting_time_in_min;
                        m_LinkCDCount[link] += 1;
                        m_LinkCACount[next_link_seq_no] += 1;

                        ////test rerouting
                        //if (bRealTimeInformationActivated)
                        //    UpdateRTPath(p_agent);
                    }

                    //move
                    p_agent->m_current_link_seq_no += 1;
                    m_LinkOutFlowCapacity[link][time_in_sec] -= 1;//deduct the capacity
                }
            }
        } // conditions
    }  // departure time events

    simu_log_file.close();

}

// FILE* g_pFileOutputLog = nullptr;
