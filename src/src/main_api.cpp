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



__int64 g_GetCellID(double x, double y, double grid_resolution)
{
    __int64 xi;
    xi = floor(x / grid_resolution);

    __int64 yi;
    yi = floor(y / grid_resolution);

    __int64 x_code, y_code, code;
    x_code = fabs(xi) * grid_resolution * 1000000000000;
    y_code = fabs(yi) * grid_resolution * 100000;
    code = x_code + y_code;
    return code;
};

string g_GetCellCode(double x, double y, double grid_resolution, double left, double top)
{
    std::string s("ABCDEFGHIJKLMNOPQRSTUVWXYZ");
    std::string str_letter;
    std::string code;

    __int64 xi;
    xi = floor(x / grid_resolution) - floor(left / grid_resolution);

    __int64 yi;
    yi = ceil(top / grid_resolution) - floor(y / grid_resolution);

    int digit = (int)(xi / 26);
    if (digit >= 1)
        str_letter = s.at(digit % s.size());

    int reminder = xi - digit * 26;
    str_letter += s.at(reminder % s.size());

    std::string num_str = std::to_string(yi);

    code = str_letter + num_str;

    return code;

}

#include "DTA.h"
#include "routing.h"
// some basic parameters setting




std::vector<NetworkForSP*> g_NetworkForSP_vector;
std::vector<NetworkForSP*> g_NetworkForRTSP_vector;
NetworkForSP g_RoutingNetwork;
std::map<int, DTAVehListPerTimeInterval> g_AgentTDListMap;
vector<CAgent_Simu*> g_agent_simu_vector;
std::vector<CNode> g_node_vector;
std::vector<CLink> g_link_vector;
std::vector<COZone> g_zone_vector;

Assignment assignment;

std::map<string, CTMC_Corridor_Info> g_tmc_corridor_vector;
std::map<string, CInfoCell> g_info_cell_map;

#include "trip_generation.h"
#include "input.h"
#include "output.h"
#include "assignment.h"
#include "ODME.h"
#include "scenario_API.h"
#include "cbi_tool.h"

std::vector<CTMC_Link> g_TMC_vector;

void g_reset_link_volume_in_master_program_without_columns(int number_of_links, int iteration_index, bool b_self_reducing_path_volume)
{
    int number_of_demand_periods = assignment.g_number_of_demand_periods;

    if(iteration_index == 0)
    {
        for (int i = 0; i < number_of_links; ++i)
        {
            for (int tau = 0; tau < number_of_demand_periods; ++tau)
            {
                // used in travel time calculation
                g_link_vector[i].flow_volume_per_period[tau] = 0;
            }
        }
    }
    else
    {
        for (int i = 0; i < number_of_links; ++i)
        {
            for (int tau = 0; tau < number_of_demand_periods; ++tau)
            {
                if (b_self_reducing_path_volume)
                {
                    // after link volumn "tally", self-deducting the path volume by 1/(k+1) (i.e. keep k/(k+1) ratio of previous flow)
                    // so that the following shortes path will be receiving 1/(k+1) flow
                    g_link_vector[i].flow_volume_per_period[tau] = g_link_vector[i].flow_volume_per_period[tau] * (float(iteration_index) / float(iteration_index + 1));
                }
            }
        }
    }
}

//***
// major function 1:  allocate memory and initialize the data
// void AllocateMemory(int number_of_nodes)
//
//major function 2: // time-dependent label correcting algorithm with double queue implementation
//int optimal_label_correcting(int origin_node, int destination_node, int departure_time, int g_debugging_flag, FILE* g_pFileDebugLog, Assignment& assignment, int time_period_no = 0, int agent_type = 1, float VOT = 10)

//	//major function: update the cost for each node at each SP tree, using a stack from the origin structure
//int tree_cost_updating(int origin_node, int departure_time, int g_debugging_flag, FILE* g_pFileDebugLog, Assignment& assignment, int time_period_no = 0, int agent_type = 1)

//***

// The one and only application object

int g_number_of_CPU_threads()
{
    int number_of_threads = omp_get_max_threads();
    int max_number_of_threads = 4000;

    if (number_of_threads > max_number_of_threads)
        number_of_threads = max_number_of_threads;

    return number_of_threads;
}

void g_assign_computing_tasks_to_memory_blocks(Assignment& assignment)
{
    //fprintf(g_pFileDebugLog, "-------g_assign_computing_tasks_to_memory_blocks-------\n");
    // step 2: assign node to thread
    dtalog.output() << "Step 2: Assigning computing tasks to memory blocks..." << endl;

    NetworkForSP* PointerMatrx[MAX_AGNETTYPES][MAX_TIMEPERIODS][MAX_MEMORY_BLOCKS] = {NULL};
    NetworkForSP* RTPointerMatrx[MAX_AGNETTYPES][MAX_TIMEPERIODS][MAX_MEMORY_BLOCKS] = { NULL };

    int computing_zone_count = 0;

    for (int at = 0; at < assignment.g_AgentTypeVector.size(); ++at)
    {
        for (int tau = 0; tau < assignment.g_DemandPeriodVector.size(); ++tau)
        {
            //assign all nodes to the corresponding thread
            for (int z = 0; z < g_zone_vector.size(); ++z)
            {

                if (z < assignment.g_number_of_memory_blocks)
                {
                    NetworkForSP* p_NetworkForSP = new NetworkForSP();

                    p_NetworkForSP->m_origin_node_vector.push_back(g_zone_vector[z].node_seq_no);
                    p_NetworkForSP->m_origin_zone_seq_no_vector.push_back(z);

                    p_NetworkForSP->m_agent_type_no = at;
                    p_NetworkForSP->m_tau = tau;
                    computing_zone_count++;

                    p_NetworkForSP->AllocateMemory(assignment.g_number_of_nodes, assignment.g_number_of_links);

                    PointerMatrx[at][tau][z] = p_NetworkForSP;

                    g_NetworkForSP_vector.push_back(p_NetworkForSP);
                }
                else  // zone seq no is greater than g_number_of_memory_blocks
                {
                    if (assignment.g_origin_demand_array[z] > 0.001 ||
             assignment.zone_seq_no_2_info_mapping.find(z) != assignment.zone_seq_no_2_info_mapping.end() 
                        )
                     {
                         //get the corresponding memory block seq no
                        // take residual of memory block size to map a zone no to a memory block no.
                         int memory_block_no = z % assignment.g_number_of_memory_blocks;
                         NetworkForSP* p_NetworkForSP = PointerMatrx[at][tau][memory_block_no];
                         p_NetworkForSP->m_origin_node_vector.push_back(g_zone_vector[z].node_seq_no);
                         p_NetworkForSP->m_origin_zone_seq_no_vector.push_back(z);
                         computing_zone_count++;
                     }
                }
            }
        }
        
    }


    dtalog.output() << "There are " << g_NetworkForSP_vector.size() << " SP networks in memory." << endl;
    dtalog.output() << "There are " << g_NetworkForRTSP_vector.size() << " RTSP networks in memory." << endl;
    dtalog.output() << "There are " << computing_zone_count << " agent type*zones to be computed in CPU." << endl;

}

void g_deallocate_computing_tasks_from_memory_blocks()
{
    //fprintf(g_pFileDebugLog, "-------g_assign_computing_tasks_to_memory_blocks-------\n");
    // step 2: assign node to thread
    for (int n = 0; n < g_NetworkForSP_vector.size(); ++n)
    {
        NetworkForSP* p_NetworkForSP = g_NetworkForSP_vector[n];
        delete p_NetworkForSP;
    }
}

//void g_reset_link_volume_for_all_processors()
//{
//#pragma omp parallel for
//    for (int ProcessID = 0; ProcessID < g_NetworkForSP_vector.size(); ++ProcessID)
//    {
//        NetworkForSP* pNetwork = g_NetworkForSP_vector[ProcessID];
//        //Initialization for all non-origin nodes
//        int number_of_links = assignment.g_number_of_links;
//        for (int i = 0; i < number_of_links; ++i)
//            pNetwork->m_link_flow_volume_array[i] = 0;
//    }
//}


void g_reset_link_volume_for_all_processors()
{
	int number_of_memory_blocks = min((int)g_NetworkForSP_vector.size(), assignment.g_number_of_memory_blocks);

#pragma omp parallel for
	for (int ProcessID = 0; ProcessID < number_of_memory_blocks; ++ProcessID)
	{
		for (int blk = 0; blk < assignment.g_AgentTypeVector.size()*assignment.g_DemandPeriodVector.size(); ++blk)
		{
			NetworkForSP* pNetwork = g_NetworkForSP_vector[blk*assignment.g_number_of_memory_blocks + ProcessID];
			//Initialization for all non-origin nodes
			int number_of_links = assignment.g_number_of_links;
			for (int i = 0; i < number_of_links; ++i)
				pNetwork->m_link_flow_volume_array[i] = 0;
		}

	}
}


void g_fetch_link_volume_for_all_processors()
{
    for (int ProcessID = 0; ProcessID < g_NetworkForSP_vector.size(); ++ProcessID)
    {
        NetworkForSP* pNetwork = g_NetworkForSP_vector[ProcessID];

        for (int i = 0; i < g_link_vector.size(); ++i)
            g_link_vector[i].flow_volume_per_period[pNetwork->m_tau] += pNetwork->m_link_flow_volume_array[i];
    }
    // step 1: travel time based on VDF
}

//major function: update the cost for each node at each SP tree, using a stack from the origin structure
double NetworkForSP::backtrace_shortest_path_tree(Assignment& assignment, int iteration_number_outterloop, int o_node_index)
{
	double total_origin_least_cost = 0;
    int origin_node = m_origin_node_vector[o_node_index]; // assigned no
    int m_origin_zone_seq_no = m_origin_zone_seq_no_vector[o_node_index]; // assigned no

    //if (assignment.g_number_of_nodes >= 100000 && m_origin_zone_seq_no % 100 == 0)
    //{
    //	g_fout << "backtracing for zone " << m_origin_zone_seq_no << endl;
    //}

    int departure_time = m_tau;
    int agent_type = m_agent_type_no;

    if (g_node_vector[origin_node].m_outgoing_link_seq_no_vector.size() == 0)
        return 0;

    // given,  m_node_label_cost[i]; is the gradient cost , for this tree's, from its origin to the destination node i'.

    //	fprintf(g_pFileDebugLog, "------------START: origin:  %d  ; Departure time: %d  ; demand type:  %d  --------------\n", origin_node + 1, departure_time, agent_type);
    float k_path_prob = 1;
    k_path_prob = float(1) / float(iteration_number_outterloop + 1);  //XZ: use default value as MSA, this is equivalent to 1/(n+1)
    // MSA to distribute the continuous flow
    // to do, this is for each nth tree.

    //change of path flow is a function of cost gap (the updated node cost for the path generated at the previous iteration -m_node_label_cost[i] at the current iteration)
    // current path flow - change of path flow,
    // new propability for flow staying at this path
    // for current shortest path, collect all the switched path from the other shortest paths for this ODK pair.
    // check demand flow balance constraints

    int num = 0;
    int number_of_nodes = assignment.g_number_of_nodes;
    int number_of_links = assignment.g_number_of_links;
    int l_node_size = 0;  // initialize local node size index
    int l_link_size = 0;
    int node_sum = 0;

    float path_travel_time = 0;
    float path_distance = 0;

    int current_node_seq_no = -1;  // destination node
    int current_link_seq_no = -1;
    int destination_zone_seq_no;
    double ODvolume, volume;
    CColumnVector* pColumnVector;

    for (int i = 0; i < number_of_nodes; ++i)
    {
        if (g_node_vector[i].zone_id >= 1)
        {
            if (i == origin_node) // no within zone demand
                continue;
            //fprintf(g_pFileDebugLog, "--------origin  %d ; destination node: %d ; (zone: %d) -------\n", origin_node + 1, i+1, g_node_vector[i].zone_id);
            //fprintf(g_pFileDebugLog, "--------iteration number outterloop  %d ;  -------\n", iteration_number_outterloop);
            destination_zone_seq_no = assignment.g_zoneid_to_zone_seq_no_mapping[g_node_vector[i].zone_id];

            pColumnVector = &(assignment.g_column_pool[m_origin_zone_seq_no][destination_zone_seq_no][agent_type][m_tau]);

            if (pColumnVector->bfixed_route) // with routing policy, no need to run MSA for adding new columns
                continue;

            ODvolume = pColumnVector->od_volume;
            volume = ODvolume * k_path_prob;
            // this is contributed path volume from OD flow (O, D, k, per time period

            if (ODvolume > 0.000001 ||
                assignment.zone_seq_no_2_activity_mapping.find(m_origin_zone_seq_no) != assignment.zone_seq_no_2_activity_mapping.end() ||
                assignment.zone_seq_no_2_activity_mapping.find(destination_zone_seq_no) != assignment.zone_seq_no_2_activity_mapping.end() ||
                (
                assignment.zone_seq_no_2_info_mapping.find(m_origin_zone_seq_no)!= assignment.zone_seq_no_2_info_mapping.end()
                && assignment.g_AgentTypeVector[agent_type].real_time_information >= 1)
                )
            {
                l_node_size = 0;  // initialize local node size index
                l_link_size = 0;
                node_sum = 0;

                path_travel_time = 0;
                path_distance = 0;

                current_node_seq_no = i;  // destination node
                current_link_seq_no = -1;

				total_origin_least_cost += ODvolume * m_node_label_cost[current_node_seq_no];
                // backtrace the sp tree from the destination to the root (at origin)
                while (current_node_seq_no >= 0 && current_node_seq_no < number_of_nodes)
                {
                    temp_path_node_vector[l_node_size++] = current_node_seq_no;

                    if (l_node_size >= temp_path_node_vector_size)
                    {
                        dtalog.output() << "Error: l_node_size >= temp_path_node_vector_size" << endl;
                        g_program_stop();
                    }

                    // this is valid node
                    if (current_node_seq_no >= 0 && current_node_seq_no < number_of_nodes)
                    {
                        current_link_seq_no = m_link_predecessor[current_node_seq_no];
                        node_sum += current_node_seq_no * current_link_seq_no;

                        // fetch m_link_predecessor to mark the link volume
                        if (current_link_seq_no >= 0 && current_link_seq_no < number_of_links)
                        {
                            temp_path_link_vector[l_link_size++] = current_link_seq_no;

                            // pure link based computing mode
                            if (assignment.assignment_mode == 0)
                            {
                                // this is critical for parallel computing as we can write the volume to data
                                m_link_flow_volume_array[current_link_seq_no]+= volume;
                            }

                            //path_travel_time += g_link_vector[current_link_seq_no].travel_time_per_period[tau];
                            //path_distance += g_link_vector[current_link_seq_no].length;
                        }
                    }
                    current_node_seq_no = m_node_predecessor[current_node_seq_no];  // update node seq no
                }
                //fprintf(g_pFileDebugLog, "\n");

                // we obtain the cost, time, distance from the last tree-k
                if(assignment.assignment_mode >=1) // column based mode
                {

                    if (iteration_number_outterloop == 1)
                    {
                        int i_debug = 2;
                    }

                    // we cannot find a path with the same node sum, so we need to add this path into the map,
                    if (pColumnVector->path_node_sequence_map.find(node_sum) == assignment.g_column_pool[m_origin_zone_seq_no][destination_zone_seq_no][agent_type][m_tau].path_node_sequence_map.end())
                    {
                        // add this unique path
                        int path_count = pColumnVector->path_node_sequence_map.size();
                        pColumnVector->path_node_sequence_map[node_sum].path_seq_no = path_count;
                        pColumnVector->path_node_sequence_map[node_sum].path_volume = 0;
                        //assignment.g_column_pool[m_origin_zone_seq_no][destination_zone_seq_no][agent_type][tau].time = m_label_time_array[i];
                        //assignment.g_column_pool[m_origin_zone_seq_no][destination_zone_seq_no][agent_type][tau].path_node_sequence_map[node_sum].path_distance = m_label_distance_array[i];
                        pColumnVector->path_node_sequence_map[node_sum].path_toll = m_node_label_cost[i];

                        pColumnVector->path_node_sequence_map[node_sum].AllocateVector(
                            l_node_size,
                            temp_path_node_vector,
                            l_link_size,
                            temp_path_link_vector,true);
                    }

                    pColumnVector->path_node_sequence_map[node_sum].path_volume += volume;
                }
            }
        }
    }
	return total_origin_least_cost;
}

void  CLink::calculate_dynamic_VDFunction(int inner_iteration_number, bool congestion_bottleneck_sensitivity_analysis_mode, int VDF_type_no)
{
    RT_travel_time = 0; // reset RT_travel time for each end of simulation iteration 
    for (int tau = 0; tau < assignment.g_number_of_demand_periods; ++tau)
    {
        float starting_time_slot_no = assignment.g_DemandPeriodVector[tau].starting_time_slot_no;
        float ending_time_slot_no = assignment.g_DemandPeriodVector[tau].ending_time_slot_no;

        float volume_to_be_assigned = flow_volume_per_period[tau];
        if (congestion_bottleneck_sensitivity_analysis_mode == true)
        {
            volume_to_be_assigned = VDF_period[tau].sav;
        }


        travel_time_per_period[tau] = VDF_period[tau].calculate_travel_time_based_on_VDF(
            volume_to_be_assigned,
            VDF_type_no, 
            this->number_of_lanes,
            starting_time_slot_no, ending_time_slot_no,
            assignment.g_DemandPeriodVector[tau].t2_peak_in_hour,
            this->free_speed, this->kc, this->s3_m, this->vc,
            this->est_speed, this->est_volume_per_hour_per_lane);

        VDF_period[tau].travel_time_per_iteration_map[inner_iteration_number] = travel_time_per_period[tau];
    }
}

double network_assignment(int assignment_mode, int iteration_number, int column_updating_iterations, int ODME_iterations, int number_of_memory_blocks)
{
    int signal_updating_iterations = 0;

    // k iterations for column generation
    assignment.g_number_of_column_generation_iterations = iteration_number;
    assignment.g_number_of_column_updating_iterations = column_updating_iterations;
    // 0: link UE: 1: path UE, 2: Path SO, 3: path resource constraints
    assignment.assignment_mode = assignment_mode;

	assignment.g_number_of_memory_blocks = min( max(1, number_of_memory_blocks), MAX_MEMORY_BLOCKS);
	
    if (assignment.assignment_mode == 0)
        column_updating_iterations = 0;

    // step 1: read input data of network / demand tables / Toll
    g_read_input_data(assignment);


    if (assignment.assignment_mode == 11)  // cbi mode
    {
        bool ReadingDataReady = assignment.map_tmc_reading();
        g_output_tmc_file(ReadingDataReady);
        return 1.0;
    }

    if (assignment.assignment_mode == 12)  // cbsa mode
    {
        int inner_iteration_number = 0;

        for (int i = 0; i < g_link_vector.size(); ++i)
        {

            // will use sav to compute travel time, time dependent speed and volume based on Q_VDF
             g_link_vector[i].VDF_type_no = 1;
             
             for (int time_index = 0; time_index < MAX_TIMEINTERVAL_PerDay; time_index++)
             {
                 g_link_vector[i].est_speed[time_index] = -1;
                 g_link_vector[i].est_volume_per_hour_per_lane[time_index] = -1;
             }
             g_link_vector[i].calculate_dynamic_VDFunction(inner_iteration_number, true, g_link_vector[i].VDF_type_no);
        }

        g_output_dynamic_QVDF_link_performance();
        return 1.0;
    }
    

    //g_reload_timing_arc_data(assignment); // no need to load timing data from timing.csv
    g_load_scenario_data(assignment);
    g_ReadDemandFileBasedOnDemandFileList(assignment);

    //step 2: allocate memory and assign computing tasks
    g_assign_computing_tasks_to_memory_blocks(assignment); // static cost based label correcting

    // definte timestamps
    clock_t start_t, end_t, total_t;
    clock_t start_simu, end_simu, total_simu;
    start_t = clock();
    clock_t iteration_t, cumulative_lc, cumulative_cp, cumulative_lu;

    //step 3: column generation stage: find shortest path and update path cost of tree using TD travel time
    dtalog.output() << endl;
    dtalog.output() << "Step 3: Column Generation for Traffic Assignment..." << endl;
    dtalog.output() << "Total Column Generation iteration: " << assignment.g_number_of_column_generation_iterations << endl;
    for (int iteration_number = 0; iteration_number < assignment.g_number_of_column_generation_iterations; iteration_number++)
    {
        dtalog.output() << endl;
        dtalog.output() << "Current iteration number:" << iteration_number << endl;
        end_t = clock();
        iteration_t = end_t - start_t;
        dtalog.output() << "Current CPU time: " << iteration_t / 1000.0 << " s" << endl;

        // step 3.1 update travel time and resource consumption
        clock_t start_t_lu = clock();

		double total_system_travel_time = 0;
		double total_least_system_travel_time = 0;
        // initialization at beginning of shortest path
		total_system_travel_time = update_link_travel_time_and_cost(iteration_number);

        if (assignment.assignment_mode == 0)
        {
            //fw
            g_reset_link_volume_in_master_program_without_columns(g_link_vector.size(), iteration_number, true);
            g_reset_link_volume_for_all_processors();
        }
        else
        {
            // we can have a recursive formulat to reupdate the current link volume by a factor of k/(k+1),
            //  and use the newly generated path flow to add the additional 1/(k+1)
            g_reset_and_update_link_volume_based_on_columns(g_link_vector.size(), iteration_number, true);
        }

        if (dtalog.debug_level() >= 3)
        {
            dtalog.output() << "Results:" << endl;
            for (int i = 0; i < g_link_vector.size(); ++i) {
                dtalog.output() << "link: " << g_node_vector[g_link_vector[i].from_node_seq_no].node_id << "-->"
                                << g_node_vector[g_link_vector[i].to_node_seq_no].node_id << ", "
                                << "flow count:" << g_link_vector[i].flow_volume_per_period[0] << endl;
            }
        }

        end_t = clock();
        iteration_t = end_t - start_t_lu;
        // g_fout << "Link update with CPU time " << iteration_t / 1000.0 << " s; " << (end_t - start_t) / 1000.0 << " s" << endl;

        //****************************************//
        //step 3.2 computng block for continuous variables;

        clock_t start_t_lc = clock();
        clock_t start_t_cp = clock();

        cumulative_lc = 0;
        cumulative_cp = 0;
        cumulative_lu = 0;

		int number_of_memory_blocks = min((int)g_NetworkForSP_vector.size(), assignment.g_number_of_memory_blocks);

#pragma omp parallel for  // step 3: C++ open mp automatically create n threads., each thread has its own computing thread on a cpu core
        //for (int ProcessID = 0; ProcessID < g_NetworkForSP_vector.size(); ++ProcessID)
        //{
        //    int agent_type_no = g_NetworkForSP_vector[ProcessID]->m_agent_type_no;

        //    for (int o_node_index = 0; o_node_index < g_NetworkForSP_vector[ProcessID]->m_origin_node_vector.size(); ++o_node_index)
        //    {
        //        start_t_lc = clock();
        //        g_NetworkForSP_vector[ProcessID]->optimal_label_correcting(ProcessID, &assignment, iteration_number, o_node_index);
        //        end_t = clock();
        //        cumulative_lc += end_t - start_t_lc;

        //        start_t_cp = clock();
        //        g_NetworkForSP_vector[ProcessID]->backtrace_shortest_path_tree(assignment, iteration_number, o_node_index);
        //        end_t = clock();
        //        cumulative_cp += end_t - start_t_cp;
        //    }
        //    // perform one to all shortest path tree calculation
        //}

		for (int ProcessID = 0; ProcessID < number_of_memory_blocks; ++ProcessID)
		{
			for (int blk = 0; blk < assignment.g_AgentTypeVector.size()*assignment.g_DemandPeriodVector.size(); ++blk)
			{
                int network_copy_no = blk* assignment.g_number_of_memory_blocks + ProcessID;
                if (network_copy_no >= g_NetworkForSP_vector.size())
                    continue;

				NetworkForSP* pNetwork = g_NetworkForSP_vector[network_copy_no];

				for (int o_node_index = 0; o_node_index < pNetwork->m_origin_node_vector.size(); ++o_node_index)
				{
					start_t_lc = clock();
					pNetwork->optimal_label_correcting(ProcessID, &assignment, iteration_number, o_node_index);


					end_t = clock();
					cumulative_lc += end_t - start_t_lc;

					start_t_cp = clock();
					double total_origin_least_travel_time = pNetwork->backtrace_shortest_path_tree(assignment, iteration_number, o_node_index);


#pragma omp critical
					{
						total_least_system_travel_time += total_origin_least_travel_time;
					}
					end_t = clock();
					cumulative_cp += end_t - start_t_cp;
				}
			}

		}

        // link based computing mode, we have to collect link volume from all processors.
        if (assignment.assignment_mode == 0)
            g_fetch_link_volume_for_all_processors();

        // g_fout << "LC with CPU time " << cumulative_lc / 1000.0 << " s; " << endl;
        // g_fout << "column generation with CPU time " << cumulative_cp / 1000.0 << " s; " << endl;

        //****************************************//

        // last iteraion before performing signal timing updating
		double relative_gap = (total_system_travel_time - total_least_system_travel_time) / max(0.00001, total_least_system_travel_time);
		dtalog.output() << "iteration: " << iteration_number << ",systemTT: " << total_system_travel_time << ", least system TT:" <<
			total_least_system_travel_time << ",gap=" << relative_gap << endl;
 }
    dtalog.output() << endl;

    // step 1.8: column updating stage: for given column pool, update volume assigned for each column
    dtalog.output() << "Step 4: Column Pool Updating" << endl;
    dtalog.output() << "Total Column Pool Updating iteration: " << column_updating_iterations << endl;
    start_t = clock();
    g_column_pool_optimization(assignment, column_updating_iterations);
    g_column_pool_route_scheduling(assignment, column_updating_iterations);
//    g_column_pool_activity_scheduling(assignment, column_updating_iterations);  // VMS information updat

    dtalog.output() << endl;

    // post route assignment aggregation
    if (assignment.assignment_mode != 0)
    {
        // we can have a recursive formulat to reupdate the current link volume by a factor of k/(k+1),
        // and use the newly generated path flow to add the additional 1/(k+1)
        g_reset_and_update_link_volume_based_on_columns(g_link_vector.size(), iteration_number, false);
    }
    else
        g_reset_link_volume_in_master_program_without_columns(g_link_vector.size(), iteration_number, false);

    // initialization at the first iteration of shortest path
    update_link_travel_time_and_cost(iteration_number);

    if (assignment.assignment_mode == 2)
    {
        start_simu = clock();

        dtalog.output() << "Step 5: Simulation for traffic assignment.." << endl;
        assignment.STTrafficSimulation();
        end_simu = clock();
        total_simu = end_simu - end_simu;

        dtalog.output() << "CPU Running Time for traffic simulation: " << total_simu / 1000.0 << " s" << endl;
        dtalog.output() << endl;
    }

    if (assignment.assignment_mode == 3)
    {
        dtalog.output() << "Step 6: O-D estimation for traffic assignment.." << endl;
        assignment.Demand_ODME(ODME_iterations);
        dtalog.output() << endl;
    }

    end_t = clock();
    total_t = (end_t - start_t);
    dtalog.output() << "Done!" << endl;

    dtalog.output() << "CPU Running Time for the entire computing progcess: " << total_t / 1000.0 << " s" << endl;

    start_t = clock();

    //step 5: output simulation results of the new demand
    g_ReadOutputFileConfiguration(assignment);
    g_output_accessibility_result(assignment);
    g_output_assignment_result(assignment);
    g_output_simulation_result(assignment);

    end_t = clock();
    total_t = (end_t - start_t);
    dtalog.output() << "Output for assignment with " << assignment.g_number_of_column_generation_iterations << " iterations. Traffic assignment completes!" << endl;
    dtalog.output() << "CPU Running Time for outputting simulation results: " << total_t / 1000.0 << " s" << endl;

    dtalog.output() << "free memory.." << endl;
    g_node_vector.clear();

    for (int i = 0; i < g_link_vector.size(); ++i)
        g_link_vector[i].free_memory();
    g_link_vector.clear();

    dtalog.output() << "done." << endl;

    return 1;
}
bool Assignment::RTSP_RealTimeShortestPathFinding(int time_slot_no, int simu_time_interval_t)
{
clock_t start_t_lc = clock();
clock_t start_t_cp = clock();

bool bRealTimeInformationActicvated = false;
for (int i = 0; i < g_link_vector.size(); ++i)
{
    CLink* pLink = &(g_link_vector[i]);
    int slot_no = time_slot_no;


}
if (bRealTimeInformationActicvated == false)  // as long as there are incident activated
return false;

for (int i = 0; i < g_link_vector.size(); ++i)
{
    CLink* pLink = &(g_link_vector[i]);
    int slot_no = time_slot_no;
    
    if (pLink->link_type >= 0)  // do not need to be the cap reduced link
           {

                    double total_waiting_time_in_min = 0;

                    for (auto it = pLink->ExitQueue.begin(); it != pLink->ExitQueue.end(); ++it)
                    {
                        int agent_id = (*it);
                        CAgent_Simu* p_agent = g_agent_simu_vector[agent_id];


                        int current_link_seq_no = p_agent->path_link_seq_no_vector[p_agent->m_current_link_seq_no];
                        int arrival_time_in_t = p_agent->m_Veh_LinkArrivalTime_in_simu_interval[p_agent->m_current_link_seq_no];
                        int waiting_time_in_t = simu_time_interval_t - arrival_time_in_t;
                        total_waiting_time_in_min += waiting_time_in_t * number_of_seconds_per_interval / 60; // in min

                    }
                    pLink->RT_travel_time = total_waiting_time_in_min / max(1, pLink->ExitQueue.size());  // average travel time

                    if (pLink->dynamic_link_closure_map.find(slot_no) != pLink->dynamic_link_closure_map.end())
                    {
                        pLink->RT_travel_time = 9999;
                    }

                    pLink->RT_speed_vector[time_slot_no] = pLink->length / max(0.00001, pLink->RT_travel_time / 60.0);  // mph                }



                    pLink->RT_travel_time_vector[time_slot_no] = pLink->RT_travel_time;

          }
        

}

 

#pragma omp parallel for  // step 3: C++ open mp automatically create n threads., each thread has its own computing thread on a cpu core
        for (int blk = 0; blk < g_NetworkForRTSP_vector.size(); ++blk)
        {

            NetworkForSP* pNetwork = g_NetworkForRTSP_vector[blk];
            int iteration_number = 0;

            for (int o_node_index = 0; o_node_index < pNetwork->m_origin_node_vector.size(); ++o_node_index)
            {
                start_t_lc = clock();
                pNetwork->optimal_label_correcting(blk, &assignment, iteration_number, o_node_index);


                start_t_cp = clock();
                double total_origin_least_travel_time = pNetwork->backtrace_shortest_path_tree(assignment, iteration_number, o_node_index);


            }
        

        }

        return true;

}

void Assignment::UpdateRTPath(CAgent_Simu* pAgent)
{
    // updating shorest path for vehicles passing through information node

    std::map<int, CColumnPath>::iterator it, it_begin, it_end;

    CColumnVector* p_column_pool;

        int at = pAgent->at;
        int dest = pAgent->dest;
        int tau = pAgent->tau;

        int current_link_seq_no = pAgent->path_link_seq_no_vector[pAgent->m_current_link_seq_no];
        CLink* pCurrentLink = &(g_link_vector[current_link_seq_no]);
        
        // step 1: check if the downstream ndoe has information zone 
        if (node_seq_no_2_info_zone_id_mapping.find(pCurrentLink->to_node_seq_no) == node_seq_no_2_info_zone_id_mapping.end())
        {
            return;
        }


        int orig = g_zoneid_to_zone_seq_no_mapping[ node_seq_no_2_info_zone_id_mapping[pCurrentLink->to_node_seq_no]];

        //step 2: fetch the related column pool from the information node/zone
        p_column_pool = &(g_column_pool[orig][dest][at][tau]);
//             scan through the map with different node sum for different continuous paths
            it_begin = p_column_pool->path_node_sequence_map.begin();
            it_end = p_column_pool->path_node_sequence_map.end();

            
            //step 3: current path has been impacted or now due to the capacity reduction

                for (it = it_begin; it != it_end; ++it)
                {
                    bool b_original_path_impacted_flag = false;
                    int link_sum_0 = 0;
                    int link_sum_rerouting = 0;
                    for (int ls = pAgent->m_current_link_seq_no+1; ls < pAgent->path_link_seq_no_vector.size(); ls++)
                    {
                        link_sum_0 += pAgent->path_link_seq_no_vector[ls];

                        if (g_link_vector[pAgent->path_link_seq_no_vector[ls]].capacity_reduction_map.find(tau) != g_link_vector[pAgent->path_link_seq_no_vector[ls]].capacity_reduction_map.end())
                            b_original_path_impacted_flag = true;
                    }

                    for (int nl = 0; nl < it->second.m_link_size-1; ++nl)  // arc a // exclude virtual link at the end;
                    {
                        link_sum_rerouting += it->second.path_link_vector[nl];
                    }

                    // detection
                    if (b_original_path_impacted_flag == false)  // this path is not impacted
                        continue; 

                    if (link_sum_0 != link_sum_rerouting)  
                    {
                 //step 4: reruting from the downstream node of the current path
                        //dtalog.output() << "route is changed for agent . " << pAgent->agent_id << endl;

                        int debug_flag = 0;

                        if(debug_flag==1)
                            {
                            for (int ls = pAgent->m_current_link_seq_no + 1; ls < pAgent->path_link_seq_no_vector.size(); ls++)
                            {
                                link_sum_0 += pAgent->path_link_seq_no_vector[ls];
                                //dtalog.output() << "current no." << ls << ":" << pAgent->path_link_seq_no_vector[ls] << ";" << endl;
                            }

                            for (int nl = 0; nl < it->second.m_link_size - 1; ++nl)  // arc a // exclude virtual link at the end;
                            {
                                link_sum_rerouting += it->second.path_link_vector[nl];
                                dtalog.output() << "rerouting no." << nl << ":" << it->second.path_link_vector[nl] << ";" << endl;
                            }
                            }

                        p_column_pool->od_volume += 1;// this is used flag

                    if(it->second.m_link_size>=2 && it->second.diverted_vehicle_map.find(pAgent->agent_id) == it->second.diverted_vehicle_map.end())  // feasible rerouting and have not been informed by this sensor yet
                    {

                        if (pAgent->diversion_flag >= 1)   // a vehicle is diverted once
                            continue; 

                        if (pAgent->agent_id == 5802)
                        {
                            int idebug = 1;
                        }

                    // step 5: change the the current routes 
                        pAgent->path_link_seq_no_vector.resize(pAgent->m_current_link_seq_no+1);  // consider virtual link
                        pAgent->m_Veh_LinkArrivalTime_in_simu_interval.resize(pAgent->m_current_link_seq_no+1);
                        pAgent->m_Veh_LinkDepartureTime_in_simu_interval.resize(pAgent->m_current_link_seq_no+1);

                        pAgent->diversion_flag += 1;
                        //expanding
                        for (int nl = 0; nl < it->second.m_link_size; ++nl)  // arc a  // we do not exclude virtual link at the end here. as the output will exclude the virtual link in trajectory.csv
                        {
                            int link_seq_no = it->second.path_link_vector[nl];
                            pAgent->path_link_seq_no_vector.push_back(link_seq_no);
                            pAgent->m_Veh_LinkArrivalTime_in_simu_interval.push_back(-1);
                            pAgent->m_Veh_LinkDepartureTime_in_simu_interval.push_back(-1);

                        }
                        it->second.path_volume += 1;

                        it->second.diverted_vehicle_map[pAgent->agent_id] = true;

                        break;
                    }

                    }
                }

  
}



