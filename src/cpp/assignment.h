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

#include "config.h"
#include "utils.h"
#include "DTA.h"

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


void g_reset_and_update_link_volume_based_on_columns(int number_of_links, int iteration_index, bool b_self_reducing_path_volume, bool b_sensitivity_analysis_flag)
{
	// record numbers
	if (b_sensitivity_analysis_flag)
	{
		for (int i = 0; i < number_of_links; ++i)
		{
			for (int tau = 0; tau < assignment.g_number_of_demand_periods; ++tau)
			{
				g_link_vector[i].VDF_period[tau].link_volume_per_iteration_map[iteration_index] = g_link_vector[i].PCE_volume_per_period[tau] + g_link_vector[i].VDF_period[tau].preload;
				// used in travel time calculation
			}
		}
	}
	for (int i = 0; i < number_of_links; ++i)
	{
		for (int tau = 0; tau < assignment.g_number_of_demand_periods; ++tau)
		{
			// used in travel time calculation
			g_link_vector[i].PCE_volume_per_period[tau] = 0;
			g_link_vector[i].person_volume_per_period[tau] = 0;

			// reserved for BPR-X
			g_link_vector[i].queue_link_distance_VDF_perslot[tau] = 0;

			for (int at = 0; at < assignment.g_AgentTypeVector.size(); ++at)
				g_link_vector[i].person_volume_per_period_per_at[tau][at] = 0;
		}
	}

	if (iteration_index >= 0)
	{
		for (int at = 0; at < assignment.g_AgentTypeVector.size(); ++at)  //m
		{

			std::map<int, CColumnPath>::iterator it;
			int zone_size = g_zone_vector.size();
			int tau_size = assignment.g_DemandPeriodVector.size();

			float link_volume_contributed_by_path_volume;

			int link_seq_no;
			double PCE_ratio = 1;
			double OCC_ratio = 1;

			int nl;

			std::map<int, CColumnPath>::iterator it_begin;
			std::map<int, CColumnPath>::iterator it_end;

			int column_vector_size;
			CColumnVector* p_column_pool;

			for (int orig = 0; orig < zone_size; ++orig)  // o
			{
				for (int dest = 0; dest < zone_size; ++dest) //d
				{
					for (int tau = 0; tau < tau_size; ++tau)  //tau
					{
						p_column_pool = &(assignment.g_column_pool[orig][dest][at][tau]);
						if (p_column_pool->od_volume > 0)
						{

							column_vector_size = p_column_pool->path_node_sequence_map.size();

							it_begin = p_column_pool->path_node_sequence_map.begin();
							it_end = p_column_pool->path_node_sequence_map.end();
							for (it = it_begin; it != it_end; ++it)
							{
								link_volume_contributed_by_path_volume = it->second.path_volume;  // assign all OD flow to this first path

								// add path volume to link volume
								for (nl = 0; nl < it->second.m_link_size; ++nl)  // arc a
								{
									link_seq_no = it->second.path_link_vector[nl];

									// MSA updating for the existing column pools
									// if iteration_index = 0; then update no flow discount is used (for the column pool case)
									PCE_ratio = g_link_vector[link_seq_no].VDF_period[tau].pce[at];  // updated on 08/16/2021 for link dependent and agent type dependent pce factor mainly for trucks 
									OCC_ratio = g_link_vector[link_seq_no].VDF_period[tau].occ[at];  // updated on 08/16/2021 for link dependent and agent type dependent pce factor mainly for trucks 
#pragma omp critical
									{
										g_link_vector[link_seq_no].PCE_volume_per_period[tau] += link_volume_contributed_by_path_volume * PCE_ratio;
										g_link_vector[link_seq_no].person_volume_per_period[tau] += link_volume_contributed_by_path_volume * OCC_ratio;
										g_link_vector[link_seq_no].person_volume_per_period_per_at[tau][at] += link_volume_contributed_by_path_volume;  // pure volume, not consider PCE
									}
								}

								// this  self-deducting action does not agents with fixed routing policies.
								if (!p_column_pool->bfixed_route && b_self_reducing_path_volume)
								{
									//after link volumn "tally", self-deducting the path volume by 1/(k+1) (i.e. keep k/(k+1) ratio of previous flow) so that the following shortes path will be receiving 1/(k+1) flow
									it->second.path_volume = it->second.path_volume * (float(iteration_index) / float(iteration_index + 1));
								}
							}
						}
					}
				}
			}
		}
	}
}

double update_link_travel_time_and_cost(int inner_iteration_number)
{
	if (assignment.assignment_mode == 2)
	{
		//compute the time-dependent delay from simulation
		//for (int l = 0; l < g_link_vector.size(); l++)
		//{
		//	float volume = assignment.m_LinkCumulativeDepartureVector[l][assignment.g_number_of_simulation_intervals - 1];  // link flow rates
		//	float waiting_time_count = 0;

		//for (int tt = 0; tt < assignment.g_number_of_simulation_intervals; tt++)
		//{
		//	waiting_time_count += assignment.m_link_TD_waiting_time[l][tt/number_of_simu_intervals_in_min];   // tally total waiting cou
		//}

		//for (int tau = 0; tau < assignment.g_DemandPeriodVector.size(); tau++)
		//{
		//	float travel_time = g_link_vector[l].free_flow_travel_time_in_min  + waiting_time_count* number_of_seconds_per_interval / max(1, volume) / 60;
		//	g_link_vector[l].travel_time_per_period[tau] = travel_time;

		//}
	}

#pragma omp parallel for
	for (int i = 0; i < g_link_vector.size(); ++i)
	{
		// step 1: travel time based on VDF
		g_link_vector[i].calculate_dynamic_VDFunction(inner_iteration_number, false, g_link_vector[i].vdf_type);

		for (int tau = 0; tau < assignment.g_DemandPeriodVector.size(); ++tau)
		{
			for (int at = 0; at < assignment.g_AgentTypeVector.size(); ++at)
			{
				float PCE_agent_type = assignment.g_AgentTypeVector[at].PCE;

				// step 2: marginal cost for SO
				g_link_vector[i].calculate_marginal_cost_for_agent_type(tau, at, PCE_agent_type);

				//if (g_debug_level  >= 3 && assignment.assignment_mode >= 2 && assignment.g_pFileDebugLog != NULL)
				//	fprintf(assignment.g_pFileDebugLog, "Update link cost: link %d->%d: tau = %d, at = %d, travel_marginal =  %.3f\n",

				//		g_node_vector[g_link_vector[l].from_node_seq_no].node_id,
				//		g_node_vector[g_link_vector[l].to_node_seq_no].node_id,
				//		tau, at,
				//		g_link_vector[l].travel_marginal_cost_per_period[tau][at]);
			}
		}
	}

	double total_network_travel_time = 0;
	for (int i = 0; i < g_link_vector.size(); ++i)
	{
		for (int tau = 0; tau < assignment.g_DemandPeriodVector.size(); ++tau)
		{
			total_network_travel_time += g_link_vector[i].VDF_period[tau].avg_travel_time * g_link_vector[i].VDF_period[tau].link_volume;
		}

	}
	return total_network_travel_time;
}


// changes here are also for odmes, don't need to implement the changes in this function for now
double g_reset_and_update_link_volume_based_on_ODME_columns(int number_of_links, int iteration_no, double& system_gap)
{
	float total_gap = 0;
	float sub_total_gap_link_count = 0;
	float sub_total_system_gap_count = 0;
	system_gap = 0;
	float sub_total_gap_P_count = 0;
	float sub_total_gap_A_count = 0;

	double total_system_travel_cost = 0;
	double total_system_travel_time = 0;
	double total_system_demand = 0;
	double total_system_UE_gap = 0;
	// reset the link volume
	for (int i = 0; i < number_of_links; ++i)
	{
		for (int tau = 0; tau < assignment.g_number_of_demand_periods; ++tau)
		{
			// used in travel time calculation
			g_link_vector[i].PCE_volume_per_period[tau] = 0;
			g_link_vector[i].person_volume_per_period[tau] = 0;

			for (int at = 0; at < assignment.g_AgentTypeVector.size(); ++at)
				g_link_vector[i].person_volume_per_period_per_at[tau][at] = 0;
		}
	}

	// reset the estimated production and attraction
	for (int orig = 0; orig < g_zone_vector.size(); ++orig)  // o
	{
		g_zone_vector[orig].est_attraction = 0;
		g_zone_vector[orig].est_production = 0;
	}

	for (int at = 0; at < assignment.g_AgentTypeVector.size(); ++at)  //m
	{

		int zone_size = g_zone_vector.size();
		int tau_size = assignment.g_DemandPeriodVector.size();
		float PCE_ratio = assignment.g_AgentTypeVector[at].PCE;
		float OCC_ratio = assignment.g_AgentTypeVector[at].OCC;


#pragma omp parallel for
		for (int orig = 0; orig < zone_size; ++orig)  // o
		{
			std::map<int, CColumnPath>::iterator it;
			float link_volume_contributed_by_path_volume;
			int nl;

			std::map<int, CColumnPath>::iterator it_begin;
			std::map<int, CColumnPath>::iterator it_end;

			int column_vector_size;
			CColumnVector* p_column_pool;
			for (int dest = 0; dest < zone_size; ++dest) //d
			{
				for (int tau = 0; tau < tau_size; ++tau)  //tau
				{
					p_column_pool = &(assignment.g_column_pool[orig][dest][at][tau]);
					if (p_column_pool->od_volume > 0)
					{


						// continuous: type 0
						column_vector_size = p_column_pool->path_node_sequence_map.size();

						it_begin = p_column_pool->path_node_sequence_map.begin();
						it_end = p_column_pool->path_node_sequence_map.end();

						double least_cost = 999999;
						int least_cost_path_seq_no = -1;
						int least_cost_path_node_sum_index = -1;
						int path_seq_count = 0;

						double path_toll = 0;
						double path_gradient_cost = 0;
						double path_distance = 0;
						double path_travel_time = 0;
						int link_seq_no;

						double link_travel_time;
						double total_switched_out_path_volume = 0;

						double step_size = 0;
						double previous_path_volume = 0;

						least_cost = 999999;
						path_seq_count = 0;

						it_begin = p_column_pool->path_node_sequence_map.begin();
						it_end = p_column_pool->path_node_sequence_map.end();
						for (it = it_begin; it != it_end; ++it)
						{
							total_system_demand += it->second.path_volume;
							path_toll = 0;
							path_gradient_cost = 0;
							path_distance = 0;
							path_travel_time = 0;

							for (int nl = 0; nl < it->second.m_link_size; ++nl)  // arc a
							{
								link_seq_no = it->second.path_link_vector[nl];
								link_travel_time = g_link_vector[link_seq_no].travel_time_per_period[tau];
								path_travel_time += link_travel_time;

							}

							it->second.path_toll = path_toll;
							it->second.path_travel_time = path_travel_time;
							total_system_travel_time += (it->second.path_travel_time * it->second.path_volume);

							if (column_vector_size == 1)  // only one path
							{

								break;
							}

							if (path_travel_time < least_cost)
							{
								least_cost = path_travel_time;
								least_cost_path_seq_no = it->second.path_seq_no;
								least_cost_path_node_sum_index = it->first;
							}
#pragma omp critical
							{
								total_system_travel_cost += (it->second.path_travel_time * it->second.path_volume);
							}
						} // end for each path



						if (column_vector_size >= 2)
						{
							// step 2: calculate gradient cost difference for each column path
							total_switched_out_path_volume = 0;
							for (it = it_begin; it != it_end; ++it)
							{
								if (it->second.path_seq_no != least_cost_path_seq_no)  //for non-least cost path
								{
									it->second.UE_gap = it->second.path_travel_time - least_cost;
									it->second.UE_relative_gap = (it->second.path_travel_time - least_cost) / max(0.0001, least_cost);

#pragma omp critical
									{
										total_system_UE_gap += (it->second.UE_gap * it->second.path_volume);
									}




								}
							}

						}  // end for each path

						for (it = it_begin; it != it_end; ++it)  // path k
						{
							link_volume_contributed_by_path_volume = it->second.path_volume;  // assign all OD flow to this first path

#pragma omp critical
							{
								g_zone_vector[orig].est_production += it->second.path_volume;
								g_zone_vector[dest].est_attraction += it->second.path_volume;
							}
							// add path volume to link volume
							for (nl = 0; nl < it->second.m_link_size; ++nl)  // arc a
							{
								link_seq_no = it->second.path_link_vector[nl];

								// MSA updating for the existing column pools
								// if iteration_index = 0; then update no flow discount is used (for the column pool case)
#pragma omp critical
								{
									g_link_vector[link_seq_no].PCE_volume_per_period[tau] += link_volume_contributed_by_path_volume * PCE_ratio;
									g_link_vector[link_seq_no].person_volume_per_period[tau] += link_volume_contributed_by_path_volume * OCC_ratio;
									g_link_vector[link_seq_no].person_volume_per_period_per_at[tau][at] += link_volume_contributed_by_path_volume;  // pure volume, not consider PCE
								}
							}
						}
					}
				}
			}
		}
	}

	int total_link_count = 0;

	// calcualte deviation for each measurement type
	for (int i = 0; i < number_of_links; ++i)
	{
		g_link_vector[i].calculate_dynamic_VDFunction(iteration_no, false, g_link_vector[i].vdf_type);
		for (int tau = 0; tau < assignment.g_DemandPeriodVector.size(); ++tau)  //tau
		{
			if (assignment.g_DemandPeriodVector[tau].number_of_demand_files == 0)
				continue;

			if (g_link_vector[i].VDF_period[tau].obs_count >= 1)  // with data
			{
				g_link_vector[i].VDF_period[tau].est_count_dev = g_link_vector[i].PCE_volume_per_period[tau] + g_link_vector[i].VDF_period[tau].preload - g_link_vector[i].VDF_period[tau].obs_count;

				if (dtalog.debug_level() == 2)
				{
					dtalog.output() << "link " << g_node_vector[g_link_vector[i].from_node_seq_no].node_id
						<< "->" << g_node_vector[g_link_vector[i].to_node_seq_no].node_id
						<< "obs:, " << g_link_vector[i].VDF_period[tau].obs_count << "est:, " << g_link_vector[i].PCE_volume_per_period[tau]
						<< "dev:," << g_link_vector[i].VDF_period[tau].est_count_dev << endl;
				}
				if (g_link_vector[i].VDF_period[tau].upper_bound_flag == 0)
				{
					total_gap += abs(g_link_vector[i].VDF_period[tau].est_count_dev);
					sub_total_gap_link_count += fabs(g_link_vector[i].VDF_period[tau].est_count_dev / g_link_vector[i].VDF_period[tau].obs_count);
					sub_total_system_gap_count += g_link_vector[i].VDF_period[tau].est_count_dev / g_link_vector[i].VDF_period[tau].obs_count;
				}
				else
				{  // upper bound constraints 
					if (g_link_vector[i].VDF_period[tau].est_count_dev > 0)
					{
						total_gap += abs(g_link_vector[i].VDF_period[tau].est_count_dev);
						sub_total_gap_link_count += fabs(g_link_vector[i].VDF_period[tau].est_count_dev / g_link_vector[i].VDF_period[tau].obs_count);
						sub_total_system_gap_count += g_link_vector[i].VDF_period[tau].est_count_dev / g_link_vector[i].VDF_period[tau].obs_count;
					}
				}
				total_link_count += 1;
			}
		}

	}
	//for (int orig = 0; orig < g_zone_vector.size(); ++orig)  // o
	//{
	//    if (g_zone_vector[orig].obs_attraction >= 1)  // with observation
	//    {
	//        g_zone_vector[orig].est_attraction_dev = g_zone_vector[orig].est_attraction - g_zone_vector[orig].obs_attraction;

	//        if (dtalog.debug_level() == 2)
	//        {
	//            dtalog.output() << "zone " << g_zone_vector[orig].zone_id << "A: obs:" << g_zone_vector[orig].obs_attraction
	//                << ",est:," << g_zone_vector[orig].est_attraction << ",dev:," << g_zone_vector[orig].est_attraction_dev << endl;
	//        }

	//        total_gap += abs(g_zone_vector[orig].est_attraction_dev);
	//        sub_total_gap_A_count += g_zone_vector[orig].est_attraction_dev / g_zone_vector[orig].obs_attraction;
	//    }

	//    if (g_zone_vector[orig].obs_production >= 1)  // with observation
	//    {
	//        g_zone_vector[orig].est_production_dev = g_zone_vector[orig].est_production - g_zone_vector[orig].obs_production;

	//        if (dtalog.debug_level() == 2)
	//        {
	//            dtalog.output() << "zone " << g_zone_vector[orig].zone_id << "P: obs:" << g_zone_vector[orig].obs_production
	//                << ",est:," << g_zone_vector[orig].est_production << ",dev:," << g_zone_vector[orig].est_production_dev << endl;
	//        }

	//        total_gap += abs(g_zone_vector[orig].est_production_dev);
	//        sub_total_gap_P_count += g_zone_vector[orig].est_production_dev / g_zone_vector[orig].obs_production;
	//    }
	//}

	dtalog.output() << "ODME #" << iteration_no
		<< ", link MAE= " << total_gap / max(1, total_link_count)
		<< ",link_MAPE: " << (sub_total_gap_link_count) / max(1, total_link_count) * 100 <<
		"%,system_MPE: " << (sub_total_system_gap_count) / max(1, total_link_count) * 100 <<
		"%,avg_tt = " << total_system_travel_time / max(0.1, total_system_demand) << "(min) " <<
		",UE gap =" << total_system_UE_gap / max(0.00001, total_system_demand) << "(min)" <<
		" = (" << total_system_UE_gap / max(0.00001, total_system_travel_time) * 100 << " %)"
		<< endl;
	double gap = sub_total_gap_link_count / max(1, total_link_count);
	system_gap = sub_total_system_gap_count / max(1, total_link_count);

	return gap;
}

void g_update_gradient_cost_and_assigned_flow_in_column_pool(Assignment& assignment, int inner_iteration_number, bool b_sensitivity_analysis_flag)
{
	double total_system_cost_gap = 0;
	float total_relative_gap = 0;
	double total_system_travel_cost = 0;
	double total_system_travel_time = 0;
	double total_system_demand = 0;

	// we can have a recursive formulat to reupdate the current link volume by a factor of k/(k+1),
	// and use the newly generated path flow to add the additional 1/(k+1)
	g_reset_and_update_link_volume_based_on_columns(g_link_vector.size(), inner_iteration_number, false, b_sensitivity_analysis_flag);

	if (b_sensitivity_analysis_flag == true)  // check estimation counts
	{
		for (int i = 0; i < g_link_vector.size(); ++i)
		{
			for (int tau = 0; tau < assignment.g_number_of_demand_periods; ++tau)
			{

				if (g_link_vector[i].VDF_period[tau].obs_count >= 1)  // with data
				{
					g_link_vector[i].VDF_period[tau].est_count_dev = g_link_vector[i].PCE_volume_per_period[tau] + g_link_vector[i].VDF_period[tau].preload - g_link_vector[i].VDF_period[tau].obs_count;
				}
			}
		}
	}
	// step 4: based on newly calculated path volumn, update volume based travel time, and update volume based resource balance, update gradie
	update_link_travel_time_and_cost(inner_iteration_number);
	// step 0


//   assignment.summary_file << ",iteration,key,o,d,at,tau,volume,"<< endl;
	//step 1: calculate shortest path at inner iteration of column flow updating
//#pragma omp parallel for
	for (int orig = 0; orig < g_zone_vector.size(); ++orig)  // o
	{
		CColumnVector* p_column_pool;
		std::map<int, CColumnPath>::iterator it, it_begin, it_end;
		int column_vector_size;

		double least_gradient_cost = 999999;
		int least_gradient_cost_path_seq_no = -1;
		int least_gradient_cost_path_node_sum_index = -1;
		int path_seq_count = 0;

		double path_toll = 0;
		double path_gradient_cost = 0;
		double path_distance = 0;
		double path_travel_time = 0;
		int link_seq_no;

		double link_travel_time;
		double total_switched_out_path_volume = 0;

		double step_size = 0;
		double previous_path_volume = 0;

		for (int dest = 0; dest < g_zone_vector.size(); ++dest) //d
		{
			for (int at = 0; at < assignment.g_AgentTypeVector.size(); ++at)  //m
			{
				for (int tau = 0; tau < assignment.g_DemandPeriodVector.size(); ++tau)  //tau
				{
					p_column_pool = &(assignment.g_column_pool[orig][dest][at][tau]);
					if (p_column_pool->od_volume > 0)
					{

						double diff = p_column_pool->od_volume - p_column_pool->prev_od_volume;

						if (b_sensitivity_analysis_flag && inner_iteration_number >= 1)
						{
							if (diff < -0.0001 || diff > 0.0001)
							{
								int idebug = 1;
							}

							if (inner_iteration_number >= 1)
								diff = p_column_pool->od_volume - p_column_pool->od_volume_per_iteration_map[inner_iteration_number - 1];

							if (diff < -0.0001 || diff > 0.0001)
							{
								int idebug = 1;
							}

						}

						if (b_sensitivity_analysis_flag)
						{
							if (g_zone_vector[orig].zone_id == 6 && g_zone_vector[dest].zone_id == 2)
							{
								int idebug = 1;
							}
						}

						p_column_pool->prev_od_volume = p_column_pool->od_volume;

						column_vector_size = p_column_pool->path_node_sequence_map.size();


						if (b_sensitivity_analysis_flag)
						{
							p_column_pool->od_volume_per_iteration_map[inner_iteration_number] = p_column_pool->od_volume;
						}

						// scan through the map with different node sum for different paths
						/// step 1: update gradient cost for each column path

						least_gradient_cost = 999999;
						least_gradient_cost_path_seq_no = -1;
						least_gradient_cost_path_node_sum_index = -1;
						path_seq_count = 0;

						it_begin = p_column_pool->path_node_sequence_map.begin();
						it_end = p_column_pool->path_node_sequence_map.end();

						bool least_path_passing_improvement_flag = false;
						for (it = it_begin; it != it_end; ++it)
						{
							path_toll = 0;
							path_gradient_cost = 0;
							path_distance = 0;
							path_travel_time = 0;

							for (int nl = 0; nl < it->second.m_link_size; ++nl)  // arc a
							{
								link_seq_no = it->second.path_link_vector[nl];
								path_toll += g_link_vector[link_seq_no].VDF_period[tau].toll[at];
								path_distance += g_link_vector[link_seq_no].link_distance_VDF;
								link_travel_time = g_link_vector[link_seq_no].travel_time_per_period[tau];
								path_travel_time += link_travel_time;

								path_gradient_cost += g_link_vector[link_seq_no].get_generalized_first_order_gradient_cost_of_second_order_loss_for_agent_type(tau, at);
							}

							it->second.path_toll = path_toll;
							it->second.path_travel_time = path_travel_time;
							it->second.path_gradient_cost = path_gradient_cost;


							if (b_sensitivity_analysis_flag == false)
								it->second.path_time_per_iteration_map[inner_iteration_number] = path_travel_time;
							else  // SA mode
								it->second.path_time_per_iteration_SA_map[inner_iteration_number] = path_travel_time;

#pragma omp critical
							{
								total_system_travel_time += (it->second.path_travel_time * it->second.path_volume);
								total_system_demand += it->second.path_volume;


								if (column_vector_size == 1)  // only one path
								{
									total_system_travel_cost += (it->second.path_gradient_cost * it->second.path_volume);

								}
							}

							if (path_gradient_cost < least_gradient_cost)
							{
								least_gradient_cost = path_gradient_cost;
								least_gradient_cost_path_seq_no = it->second.path_seq_no;
								least_gradient_cost_path_node_sum_index = it->first;
								if (it->second.network_design_flag)
								{
									least_path_passing_improvement_flag = 1;
								}
							}
						}

						if (column_vector_size >= 2)
						{
							// step 2: calculate gradient cost difference for each column path
							total_switched_out_path_volume = 0;
							for (it = it_begin; it != it_end; ++it)
							{
								if (it->second.path_seq_no != least_gradient_cost_path_seq_no)  //for non-least cost path
								{
									it->second.path_gradient_cost_difference = it->second.path_gradient_cost - least_gradient_cost;

									//if(it->second.path_gradient_cost_difference >0.0001f)
									{
										it->second.path_gradient_cost_relative_difference = it->second.path_gradient_cost_difference / max(0.0001, least_gradient_cost);
									}

#pragma omp critical
									{
										total_system_cost_gap += (it->second.path_gradient_cost_difference * it->second.path_volume);
										total_system_travel_cost += (it->second.path_gradient_cost * it->second.path_volume);
									}

									if (b_sensitivity_analysis_flag == true)  // SA stages
									{

										//float est_count_dev = 0;
										//bool network_design_flag = false;
										//for (int nl = 0; nl < it->second.m_link_size; ++nl)  // arc a
										//{
										//    // step 3.3 link flow gradient
										//    link_seq_no = it->second.path_link_vector[nl];

										//    //if (g_link_vector[link_seq_no].tmc_corridor_name .size() > 0)
										//    //    network_design_flag = true;
										//    if (g_link_vector[link_seq_no].VDF_period[tau].obs_count >= 1)
										//    {
										//            path_gradient_cost += g_link_vector[link_seq_no].VDF_period[tau].est_count_dev;
										//            est_count_dev += g_link_vector[link_seq_no].VDF_period[tau].est_count_dev;

										//            //if (g_link_vector[link_seq_no].VDF_period[tau].network_design_flag==0 && g_link_vector[link_seq_no].VDF_period[tau].est_count_dev < 0)  // if under-report traffic 
										//            //{
										//            //    double weight_on_count = 0.0;
										//            //    it->second.path_gradient_cost_relative_difference -= weight_on_count* g_link_vector[link_seq_no].VDF_period[tau].est_count_dev;
										//            //}

										//   }
										//}

										//step_size = 0.00;

										//if (least_path_passing_improvement_flag)
										//{
//                                            if(network_design_flag == false)
										step_size = 0.05; // small changes 
								//}


//                                          step_size = 1.0 / (inner_iteration_number + 2) * p_column_pool->od_volume;
										  //if (network_design_flag)
										  //{
										  //    //                                        step_size = 1.0 / (inner_iteration_number + 2) * p_column_pool->od_volume;
										  //    assignment.summary_file << "," << inner_iteration_number
										  //        << "," << orig
										  //        << "-" << dest
										  //        << "-" << at
										  //        << "-" << tau
										  //        << "," << orig
										  //        << "," << dest
										  //        << "," << at
										  //        << "," << tau
										  //        << "," << p_column_pool->od_volume
										  //        << "," << step_size * it->second.path_gradient_cost_relative_difference
										  //        << endl;
										  //}
									}
									else
									{  // column updating step size
										step_size = 1.0 / (inner_iteration_number + 2) * p_column_pool->od_volume;
									}

									previous_path_volume = it->second.path_volume; //b

									double flow_shift = step_size * max(0.0000, it->second.path_gradient_cost_relative_difference);  //c, must be positive

									if (flow_shift > it->second.path_volume * 0.5)
									{
										flow_shift = it->second.path_volume * 0.5;
									}

									if (flow_shift >= 0.000001)
									{
										int idebug = 1;
									}

									//recall that it->second.path_gradient_cost_difference >=0
									// step 3.1: shift flow from nonshortest path to shortest path
									it->second.path_volume = max(0.0, it->second.path_volume - flow_shift);
									//d
									// 
									//we use min(step_size to ensure a path is not switching more than 1/n proportion of flow
									it->second.path_switch_volume = (previous_path_volume - it->second.path_volume);

									// d-b
									// should be nonnegative
									total_switched_out_path_volume += (previous_path_volume - it->second.path_volume);

									if (fabs(total_switched_out_path_volume) > 0.00001)
									{
										int idebug = 1;
									}
								}
							}

							//step 3.2 consider least cost path, receive all volume shifted from non-shortest path
							if (least_gradient_cost_path_seq_no != -1 && p_column_pool->path_node_sequence_map.find(least_gradient_cost_path_node_sum_index) != p_column_pool->path_node_sequence_map.end())
							{
								if (least_gradient_cost_path_node_sum_index < 100)
								{
									int i_debug = 1;
								}
								p_column_pool->path_node_sequence_map[least_gradient_cost_path_node_sum_index].path_volume += total_switched_out_path_volume;

								if (b_sensitivity_analysis_flag == false)
									p_column_pool->path_node_sequence_map[least_gradient_cost_path_node_sum_index].path_volume_per_iteration_map[inner_iteration_number] = p_column_pool->path_node_sequence_map[least_gradient_cost_path_node_sum_index].path_volume;
								else
									p_column_pool->path_node_sequence_map[least_gradient_cost_path_node_sum_index].path_volume_per_iteration_SA_map[inner_iteration_number] = p_column_pool->path_node_sequence_map[least_gradient_cost_path_node_sum_index].path_volume;
#pragma omp critical
								{

									total_system_travel_cost += (p_column_pool->path_node_sequence_map[least_gradient_cost_path_node_sum_index].path_gradient_cost *
										p_column_pool->path_node_sequence_map[least_gradient_cost_path_node_sum_index].path_volume);
								}
							}

						}
						// record path flow for all paths( including shortst path and non_shortest path)
						for (it = it_begin; it != it_end; ++it)
						{

							if (b_sensitivity_analysis_flag == false)
								it->second.path_volume_per_iteration_map[inner_iteration_number] = it->second.path_volume;
							else  //SA mode
								it->second.path_volume_per_iteration_SA_map[inner_iteration_number] = it->second.path_volume;
						}


					}
				}
			}
		}
	}

	double avg_travel_time = total_system_travel_time / max(0.001, total_system_demand);

	dtalog.output() << "column updating: iteration= " << inner_iteration_number << ", avg travel time = " <<
		avg_travel_time << "(min), optimization obj = " << total_system_cost_gap
		<< ",Relative_gap=" << total_system_cost_gap * 100.0 / max(0.00001, total_system_travel_cost) << " %" << endl;

	string stage_str;
	stage_str = "column updating";

	if (b_sensitivity_analysis_flag)
		stage_str = "sensitivity analaysis";

	assignment.summary_file2 << stage_str.c_str() << ",iteration," << inner_iteration_number <<
		",total_system_demand," << total_system_demand <<
		",avg travel time," << avg_travel_time <<
		",optimization obj," << total_system_cost_gap <<
		",relative_gap," << total_system_cost_gap * 100.0 / max(0.00001, total_system_travel_cost) << "," << endl;
}

void g_classification_in_column_pool(Assignment& assignment)
{
	int impact_OD_counts = 0;
	int impact_OD_counts_detour = 0;
	//#pragma omp parallel for
	for (int orig = 0; orig < g_zone_vector.size(); ++orig)  // o
	{
		CColumnVector* p_column_pool;
		std::map<int, CColumnPath>::iterator it, it_begin, it_end;
		int column_vector_size;

		int link_seq_no;

		for (int dest = 0; dest < g_zone_vector.size(); ++dest) //d
		{
			for (int at = 0; at < assignment.g_AgentTypeVector.size(); ++at)  //m
			{
				for (int tau = 0; tau < assignment.g_DemandPeriodVector.size(); ++tau)  //tau
				{

					p_column_pool = &(assignment.g_column_pool[orig][dest][at][tau]);

					if (p_column_pool->od_volume > 0)
					{

						if (g_zone_vector[orig].zone_id == 6 && g_zone_vector[dest].zone_id == 2)
						{
							int idebug = 1;
						}
						column_vector_size = p_column_pool->path_node_sequence_map.size();


						// scan through the map with different node sum for different paths
						/// step 1: update gradient cost for each column path

						it_begin = p_column_pool->path_node_sequence_map.begin();
						it_end = p_column_pool->path_node_sequence_map.end();

						bool least_path_passing_improvement_flag = false;
						// scan all paths in this OD pair
						int path_count = 0;
						int network_design_path_count = 0;

						for (it = it_begin; it != it_end; ++it)
						{
							for (int nl = 0; nl < it->second.m_link_size; ++nl)  // arc a
							{
								link_seq_no = it->second.path_link_vector[nl];

								if (g_link_vector[link_seq_no].VDF_period[tau].network_design_flag != 0)  // screening condition 1: passing through the network design location
								{
									it->second.network_design_flag = 1;

									// to be revised: passing through work zone, and with signal timing enhancemnets
								}

							}

							if (it->second.network_design_flag)
								network_design_path_count++;

							path_count++;
						}


						if (network_design_path_count >= 1)
						{
							if (network_design_path_count == path_count)
							{
								p_column_pool->OD_network_design_flag = 1;
								impact_OD_counts++;
							}
							else
							{
								p_column_pool->OD_network_design_flag = 2;  // more than 2 alterantive paths with respect to the newtork design location
								impact_OD_counts_detour++;
							}

						}

						if (p_column_pool->OD_network_design_flag == 2)  // 
						{
							// scan all paths in this OD pair again // mark alternative paths
							for (it = it_begin; it != it_end; ++it)
							{
								if (it->second.network_design_flag == 0)
								{
									it->second.network_design_detour_mode = 2; // detour
								}
								else
								{
									it->second.network_design_detour_mode = 1; // main passing path
								}

							}
						}

					}
				} // for each tau
			}// for each agent type mode
		} // for each d
	}
	string stage_str;
	stage_str = "classification";


//	assignment.summary_file2 << stage_str.c_str() << ",impact_OD_counts," << impact_OD_counts <<
//		",impact_OD_counts_with_detour," << impact_OD_counts_detour << endl;



}

void g_column_pool_optimization(Assignment& assignment, int column_updating_iterations, bool sensitivity_analysis_flag = false)
{
	// column_updating_iterations is internal numbers of column updating
	for (int n = 0; n < column_updating_iterations; ++n)
	{
		g_update_gradient_cost_and_assigned_flow_in_column_pool(assignment, n, sensitivity_analysis_flag);

		if (dtalog.debug_level() >= 3)
		{
			for (int i = 0; i < g_link_vector.size(); ++i)
			{
				dtalog.output() << "link: " << g_node_vector[g_link_vector[i].from_node_seq_no].node_id << "-->"
					<< g_node_vector[g_link_vector[i].to_node_seq_no].node_id << ", "
					<< "flow count:" << g_link_vector[i].PCE_volume_per_period[0] << endl;
			}
		}
	}
}


void g_column_pool_route_scheduling(Assignment& assignment, int inner_iteration_number)
{

	//step 1: calculate shortest path at inner iteration of column flow updating
#pragma omp parallel for
	for (int orig = 0; orig < g_zone_vector.size(); ++orig)  // o
	{
		CColumnVector* p_column_pool;
		std::map<int, CColumnPath>::iterator it, it_begin, it_end;
		int column_vector_size;

		int path_seq_count = 0;

		double path_toll = 0;
		double path_gradient_cost = 0;
		double path_distance = 0;
		double path_travel_time = 0;
		int link_seq_no;


		for (int dest = 0; dest < g_zone_vector.size(); ++dest) //d
		{
			for (int at = 0; at < assignment.g_AgentTypeVector.size(); ++at)  //m
			{
				for (int tau = 0; tau < assignment.g_DemandPeriodVector.size(); ++tau)  //tau
				{
					p_column_pool = &(assignment.g_column_pool[orig][dest][at][tau]);
					if (p_column_pool->od_volume > 0)
					{

						if (assignment.g_AgentTypeVector[at].real_time_information == 1)   // case of VMS 
						{

							column_vector_size = p_column_pool->path_node_sequence_map.size();

							// scan through the map with different node sum for different paths

							path_seq_count = 0;

							it_begin = p_column_pool->path_node_sequence_map.begin();
							it_end = p_column_pool->path_node_sequence_map.end();

							//test condition 1: passing through information zone
							bool b_passing_information_zone = false;
							int new_orig_zone_id = 0;

							std::vector <int> link_seq_vector;

							//test condition 2: passing through capacity impact area
							bool b_passing_capacity_impact_area = false;
							for (it = it_begin; it != it_end; ++it)  // scan each first-stage original path
							{
								if (it->second.path_volume < 0.00001)
									continue;

								for (int nl = 0; nl < it->second.m_link_size; ++nl)  // arc a
								{
									link_seq_no = it->second.path_link_vector[nl];
									CLink* p_current_link = &(g_link_vector[link_seq_no]);



									if (b_passing_information_zone == false &&
										assignment.node_seq_no_2_info_zone_id_mapping.find(p_current_link->to_node_seq_no) != assignment.node_seq_no_2_info_zone_id_mapping.end())  // this node been defined as zone
									{
										int zone_id = assignment.node_seq_no_2_info_zone_id_mapping[p_current_link->to_node_seq_no];
										int zone_no = assignment.g_zoneid_to_zone_seq_no_mapping[zone_id];

										if (assignment.zone_seq_no_2_info_mapping.find(zone_no) != assignment.zone_seq_no_2_info_mapping.end())  // as information zone
										{
											b_passing_information_zone = true;
											new_orig_zone_id = zone_id;  // zone id to zone no.

											for (int nl2 = 0; nl2 <= nl; ++nl2)  // arc a
											{  // copy the existing link sequence up to the downstream node id corresponding to the info zone
												link_seq_no = it->second.path_link_vector[nl2];
												link_seq_vector.push_back(link_seq_no);
											}
										}
									}

									if (p_current_link->capacity_reduction_map.find(tau) != p_current_link->capacity_reduction_map.end())
									{
										b_passing_capacity_impact_area = true;
									}

								}



								if (b_passing_capacity_impact_area == true && b_passing_information_zone == true)
								{

									CColumnVector* p_2_stage_column_pool;

									int info_orig = assignment.g_zoneid_to_zone_seq_no_mapping[new_orig_zone_id];

									//step 2: fetch the related column pool from the information node/zone
									p_2_stage_column_pool = &(assignment.g_column_pool[info_orig][dest][at][tau]);  // we come from info_orig but going to  the same destination with same at, and assignment period tau 
									//             scan through the map with different node sum for different continuous paths

									std::map<int, CColumnPath>::iterator it2, it_begin2, it_end2;

									it_begin2 = p_2_stage_column_pool->path_node_sequence_map.begin();
									it_end2 = p_2_stage_column_pool->path_node_sequence_map.end();


									for (it2 = it_begin2; it2 != it_end2; ++it2)  // we can still have k-path from the info zone to to final destination so we need to random select one 
									{

										for (int nl = 1; nl < it2->second.m_link_size; ++nl)  // arc a // exclude virtual link at the end;
										{
											link_seq_vector.push_back(it2->second.path_link_vector[nl]);
										}
										break; // only connect with the first available second stage path
									}

									if (it->second.path_link_vector != NULL)
									{
										// copy the updated path (stage1 + stage 2) back to the path link vector 
										delete it->second.path_link_vector;
										it->second.path_link_vector = new int[link_seq_vector.size()];
										for (int l = 0; l < link_seq_vector.size(); l++)
										{
											it->second.path_link_vector[l] = link_seq_vector[l];
										}
										it->second.m_link_size = link_seq_vector.size();

										// copy the updated path (stage1 + stage 2) back to the path node vector 
										delete it->second.path_node_vector;
										it->second.path_node_vector = new int[link_seq_vector.size() + 1];


										// first node
										it->second.path_node_vector[0] = g_link_vector[link_seq_vector[0]].from_node_seq_no;


										// remaining nodes to the end of path

										for (int l = 0; l < link_seq_vector.size(); l++)
										{
											it->second.path_node_vector[l + 1] = g_link_vector[link_seq_vector[l]].to_node_seq_no;
										}
										it->second.m_node_size = link_seq_vector.size() + 1;
									}

									p_2_stage_column_pool->od_volume += it->second.path_volume;// carry over the switching path flow to the second path volume count

									p_2_stage_column_pool->information_type = 1;
									it2->second.path_volume += it->second.path_volume;// carry over the switching path flow to the second path volume count


								}  // two conditions satisified

							}  //end of scanning for the first stage path in the column pool

						} // agent type is real time agent type

					}  // with positve OD volume


				} // tau
			}  //agent type
		} //dest
	}  // orig


	dtalog.output() << " updating";


}


void g_rt_info_column_generation(Assignment* p_assignment, float current_time_in_min, int recording_flag = 0)
{
	//dtalog.output() << "Begin the computing of " << g_NetworkForRTSP_vector.size() << " RTSP networks in CPU." << endl;

	clock_t start_t0, end_t0, total_t0;

	start_t0 = clock();
#pragma omp parallel for  // step 3: C++ open mp automatically create n threads., each thread has its own computing thread on a cpu core
	for (int blk = 0; blk < g_NetworkForRTSP_vector.size(); ++blk)
	{


		NetworkForSP* pNetwork = g_NetworkForRTSP_vector[blk];

		if (assignment.g_DemandPeriodVector[pNetwork->m_tau].starting_time_slot_no * MIN_PER_TIMESLOT > current_time_in_min)  // RT network is for a later time interval
			continue;

		pNetwork->optimal_backward_label_correcting_from_destination(blk, p_assignment, current_time_in_min, pNetwork->m_RT_dest_zone, pNetwork->m_RT_dest_node, -1, recording_flag);

	}
	end_t0 = clock();

	total_t0 = (end_t0 - start_t0);
	int second = total_t0 / 1000.0;
	int min = second / 60;
	int sec = second - min * 60;
	//dtalog.output() << "CPU Running Time for RT shortest path: " << min << " min " << sec << " sec" << endl;

}


void g_column_pool_activity_scheduling(Assignment& assignment, int inner_iteration_number)
{

	//step 1: calculate shortest path at inner iteration of column flow updating

	for (int orig = 0; orig < g_zone_vector.size(); ++orig)  // o
	{
		CColumnVector* p_column_pool;
		int path_seq_count = 0;

		double path_toll = 0;
		double path_gradient_cost = 0;
		double path_distance = 0;
		double path_travel_time = 0;

		for (int dest = 0; dest < g_zone_vector.size(); ++dest) //d
		{
			for (int at = 0; at < assignment.g_AgentTypeVector.size(); ++at)  //m
			{
				for (int tau = 0; tau < assignment.g_DemandPeriodVector.size(); ++tau)  //tau
				{
					p_column_pool = &(assignment.g_column_pool[orig][dest][at][tau]);
					if (p_column_pool->od_volume > 0)
					{
						if (p_column_pool->activity_zone_no_vector.size())   // case of activity zones
						{
							p_column_pool->path_node_sequence_map.clear();  // remove existing single OD pair based routes


							std::vector <int> link_seq_vector;
							// for each origin and detination pair in activity zone no to perform routing continuously

							for (int az = 0; az < p_column_pool->activity_zone_no_vector.size() - 1; az++) // key step: go through each activty OD pair
							{ // 0 will the origin
								// last one will destination
								int aat = p_column_pool->activity_agent_type_no_vector[az];

								CColumnVector* p_2_stage_column_pool;

								int activity_orig = p_column_pool->activity_zone_no_vector[az];
								int activity_dest = p_column_pool->activity_zone_no_vector[az + 1];

								//step 2: fetch the related column pool from the information node/zone
								p_2_stage_column_pool = &(assignment.g_column_pool[activity_orig][activity_dest][aat][tau]);  // we come from info_orig but going to  the same destination with same at, and assignment period tau 
								//             scan through the map with different node sum for different continuous paths

								std::map<int, CColumnPath>::iterator it2, it_begin2, it_end2;

								it_begin2 = p_2_stage_column_pool->path_node_sequence_map.begin();
								it_end2 = p_2_stage_column_pool->path_node_sequence_map.end();


								for (it2 = it_begin2; it2 != it_end2; ++it2)  // we can still have k-path from the info zone to to final destination so we need to random select one 
								{

									for (int nl = 1; nl < it2->second.m_link_size - 1; ++nl)  // arc a // exclude virtual link in the beginning and at the end;
									{
										link_seq_vector.push_back(it2->second.path_link_vector[nl]);
									}
									break; // only connect with the first available second stage path
								}


							}

							if (link_seq_vector.size() == 0)
							{
								int i_debug = 1;
								continue;
							}

							int node_sum = 0;
							for (int l = 0; l < link_seq_vector.size(); l++)
							{
								node_sum += link_seq_vector[l];
							}


							// add this unique path  // later we can add k activity paths
							int path_count = p_column_pool->path_node_sequence_map.size();
							p_column_pool->path_node_sequence_map[node_sum].path_seq_no = path_count;
							p_column_pool->path_node_sequence_map[node_sum].path_volume = p_column_pool->od_volume;
							p_column_pool->path_node_sequence_map[node_sum].path_toll = 0;


							p_column_pool->path_node_sequence_map[node_sum].path_link_vector = new int[link_seq_vector.size()];
							p_column_pool->path_node_sequence_map[node_sum].path_node_vector = new int[link_seq_vector.size() + 1];

							for (int l = 0; l < link_seq_vector.size(); l++)
							{
								p_column_pool->path_node_sequence_map[node_sum].path_link_vector[l] = link_seq_vector[l];
								p_column_pool->path_node_sequence_map[node_sum].path_link_STL_vector.push_back(link_seq_vector[l]);
							}
							p_column_pool->path_node_sequence_map[node_sum].m_link_size = link_seq_vector.size();

							// copy the updated path (stage1 + stage 2) back to the path node vector 

							// first node
							p_column_pool->path_node_sequence_map[node_sum].path_node_vector[0] = g_link_vector[link_seq_vector[0]].from_node_seq_no;

							// remaining nodes to the end of path

							for (int l = 0; l < link_seq_vector.size(); l++)
							{
								p_column_pool->path_node_sequence_map[node_sum].path_node_vector[l + 1] = g_link_vector[link_seq_vector[l]].to_node_seq_no;
							}
							p_column_pool->path_node_sequence_map[node_sum].m_node_size = link_seq_vector.size() + 1;
						} //end of conditions for activity chain
					}  // with positve OD volume


				} // tau
			}  //agent type
		} //dest
	}  // orig


	dtalog.output() << " updating";
}
