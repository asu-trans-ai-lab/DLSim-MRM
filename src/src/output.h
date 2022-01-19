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

void g_output_dynamic_queue_profile()  // generated from VDF, from numerical queue evolution calculation
{
    dtalog.output() << "writing link_queue_profile.csv.." << endl;

    int b_debug_detail_flag = 0;
    FILE* g_pFileLinkMOE = nullptr;

    string file_name = "link_queue_profile.csv";

    fopen_ss(&g_pFileLinkMOE, file_name.c_str(), "w");

    if (!g_pFileLinkMOE)
    {
        dtalog.output() << "File " << file_name.c_str() << " cannot be opened." << endl;
        g_program_stop();
    }
    else
    {

        // Option 2: BPR-X function
        fprintf(g_pFileLinkMOE, "link_id,tmc_corridor_name,tmc_road_sequence,tmc,link_type_name,from_node_id,to_node_id,geometry,");

        fprintf(g_pFileLinkMOE, "link_type_code,FT,AT,nlanes,link_distance_VDF,free_speed,capacity,k_critical,v_cutoff,");
        for (int tau = 0; tau < min(3, assignment.g_DemandPeriodVector.size()); tau++)
        {
            fprintf(g_pFileLinkMOE, "%s_Volume,%s_speed_BPR,%s_speed_QVDF,%s_t0,%s_t3,%s_P, %s_D,%s_DC_ratio,%s_mu,",
                assignment.g_DemandPeriodVector[tau].demand_period.c_str(),
                assignment.g_DemandPeriodVector[tau].demand_period.c_str(),
                assignment.g_DemandPeriodVector[tau].demand_period.c_str(),
                assignment.g_DemandPeriodVector[tau].demand_period.c_str(),
                assignment.g_DemandPeriodVector[tau].demand_period.c_str(),
                assignment.g_DemandPeriodVector[tau].demand_period.c_str(),
                assignment.g_DemandPeriodVector[tau].demand_period.c_str(),
                assignment.g_DemandPeriodVector[tau].demand_period.c_str());
        }

        fprintf(g_pFileLinkMOE, "");

        // hourly data
        for (int t = 6 * 60; t < 20 * 60; t += 60)
        {
            int hour = t / 60;
            int minute = t - hour * 60;

            fprintf(g_pFileLinkMOE, "vh%02d,", hour);
        }

        for (int t = 6 * 60; t < 20 * 60; t += 60)
        {
            int hour = t / 60;
            int minute = t - hour * 60;

            fprintf(g_pFileLinkMOE, "vrh%02d,", hour);
        }

        for (int t = 6 * 60; t < 20 * 60; t += 15)
        {
            int hour = t / 60;
            int minute = t - hour * 60;
            fprintf(g_pFileLinkMOE, "vr%02d:%02d,", hour, minute);
        }

        for (int t = 6 * 60; t < 20 * 60; t += 15)
        {
            int hour = t / 60;
            int minute = t - hour * 60;
            fprintf(g_pFileLinkMOE, "v%02d:%02d,", hour, minute);
        }

        for (int t = 6 * 60; t < 20 * 60; t += 5)
        {
            int hour = t / 60;
            int minute = t - hour * 60;
            fprintf(g_pFileLinkMOE, "v%02d:%02d,", hour, minute);
        }


        fprintf(g_pFileLinkMOE, "\n");

        //Initialization for all nodes
        for (int i = 0; i < g_link_vector.size(); ++i)
        {
            // virtual connectors
            if (g_link_vector[i].link_type == -1)
                continue;

            fprintf(g_pFileLinkMOE, "%s,%s,%d,%s,%s,%d,%d,\"%s\",",
                g_link_vector[i].link_id.c_str(),
                g_link_vector[i].tmc_corridor_name.c_str(),
                g_link_vector[i].tmc_road_sequence,
                g_link_vector[i].tmc_code.c_str(),

                g_link_vector[i].link_type_name.c_str(),

                g_node_vector[g_link_vector[i].from_node_seq_no].node_id,
                g_node_vector[g_link_vector[i].to_node_seq_no].node_id,
                g_link_vector[i].geometry.c_str());


            fprintf(g_pFileLinkMOE, "%s,%d,%d,%d,%f,%f,%f,%f,%f,",
                g_link_vector[i].link_type_code.c_str(),
                g_link_vector[i].FT,
                g_link_vector[i].AT,
                g_link_vector[i].number_of_lanes,
                g_link_vector[i].link_distance_VDF,
                g_link_vector[i].free_speed,
                g_link_vector[i].lane_capacity,
                g_link_vector[i].k_critical,
                g_link_vector[i].v_cutoff);

            //AM g_link_vector[i].VDF_period[0].

            int period_index = 0;

            double assignment_PMT;
            double assignment_PHT;
            double assignment_PSDT;
            double assignment_VCDT = 0;

            double assignment_VMT;
            double assignment_VHT;

            for (int tau = 0; tau < min(3, assignment.g_DemandPeriodVector.size()); tau++)
            {
            double vehicle_volume = g_link_vector[i].vehicle_flow_volume_per_period[period_index] + g_link_vector[i].vehicle_flow_volume_per_period[period_index] + g_link_vector[i].VDF_period[period_index].sa_volume;
            double person_volume = g_link_vector[i].person_flow_volume_per_period[period_index] + g_link_vector[i].vehicle_flow_volume_per_period[period_index] + g_link_vector[i].VDF_period[period_index].sa_volume;

            assignment_VMT = g_link_vector[i].link_distance_VDF * vehicle_volume;
            assignment_VHT = g_link_vector[i].travel_time_per_period[period_index]* vehicle_volume / 60.0;  // 60.0 converts min to hour

            assignment_PMT = g_link_vector[i].link_distance_VDF * person_volume;
            assignment_PHT = g_link_vector[i].travel_time_per_period[period_index] * person_volume / 60.0;  // 60.0 converts min to hour

            assignment_PSDT = (g_link_vector[i].travel_time_per_period[period_index]- g_link_vector[i].free_flow_travel_time_in_min) * person_volume / 60.0;  // 60.0 converts min to hour

            double VCTT = g_link_vector[i].link_distance_VDF / max(1, g_link_vector[i].v_cutoff) * 60;
            assignment_VCDT = max(0,g_link_vector[i].travel_time_per_period[period_index] - VCTT) * g_link_vector[i].vehicle_flow_volume_per_period[period_index] / 60.0;  // 60.0 converts min to hour


            fprintf(g_pFileLinkMOE, "%f,%f, %f,%f,%f, %f,%f,%f,%f, ",
                g_link_vector[i].vehicle_flow_volume_per_period[period_index] + g_link_vector[i].VDF_period[period_index].preload + g_link_vector[i].VDF_period[period_index].sa_volume,
                g_link_vector[i].VDF_period[period_index].avg_speed_BPR,
                g_link_vector[i].VDF_period[period_index].avg_queue_speed,

                g_link_vector[i].VDF_period[period_index].t0,
                g_link_vector[i].VDF_period[period_index].t3,
                g_link_vector[i].VDF_period[period_index].P,

                g_link_vector[i].VDF_period[period_index].lane_based_D,
                g_link_vector[i].VDF_period[period_index].DOC);
            }

            for (int t = 6 * 60; t < 20 * 60; t += 60)
            {
                float speed = g_link_vector[i].get_est_hourly_speed(t);
                fprintf(g_pFileLinkMOE, "%.3f,", speed);
            }

            for (int t = 6 * 60; t < 20 * 60; t += 60)
            {
                float speed_ratio = g_link_vector[i].get_est_hourly_speed(t) /  max(1, g_link_vector[i].v_cutoff);
                if (speed_ratio > 1)
                    speed_ratio = 1;

                fprintf(g_pFileLinkMOE, "%.3f,", speed_ratio);
            }


            for (int t = 6 * 60; t < 20 * 60; t += 15)
            {
                int time_interval = t / 5;
                float speed = g_link_vector[i].est_speed[time_interval];
                fprintf(g_pFileLinkMOE, "%.3f,", speed);
            }

            for (int t = 6 * 60; t < 20 * 60; t += 15)
            {
                int time_interval = t / 5;
                float speed_ratio = g_link_vector[i].est_speed[time_interval] / max(1, g_link_vector[i].v_cutoff);
                if (speed_ratio > 1)
                    speed_ratio = 1;

                fprintf(g_pFileLinkMOE, "%.3f,", speed_ratio);
            }

            for (int t = 6 * 60; t < 20 * 60; t += 5)
            {
                float speed = g_link_vector[i].get_est_hourly_speed(t);
                fprintf(g_pFileLinkMOE, "%.3f,", speed);
            }

            //for (int t = 6 * 60; t < 20 * 60; t += 5)
            //{
            //    int time_interval = t / 5;
            //    float waiting_ime = g_link_vector[i].est_avg_waiting_time_in_min[time_interval];
            //    fprintf(g_pFileLinkMOE, "%.2f,", waiting_ime);
            //}
            //for (int t = 6 * 60; t < 20 * 60; t += 60)
            //{
            //    float volume = g_link_vector[i].get_est_hourly_volume(t);
            //    fprintf(g_pFileLinkMOE, "%.3f,", volume);
            //}

            fprintf(g_pFileLinkMOE, "\n");

        }  // for each link l
        fclose(g_pFileLinkMOE);
    }//assignment mode 2 as simulation

}

void g_output_assignment_result(Assignment& assignment)
{
    g_OutputModelFiles(10); // label cost tree

    dtalog.output() << "writing link_performance.csv.." << endl;

    int b_debug_detail_flag = 0;
    FILE* g_pFileLinkMOE = nullptr;

    fopen_ss(&g_pFileLinkMOE, "link_performance.csv", "w");
    if (!g_pFileLinkMOE)
    {
        dtalog.output() << "File link_performance.csv cannot be opened." << endl;
        g_program_stop();
    }
    else
    {

        // Option 2: BPR-X function
        fprintf(g_pFileLinkMOE, "link_id,from_node_id,to_node_id,link_type_name,link_type_code,time_period,volume,person_volume,travel_time,speed,speed_ratio,VOC,DOC,capacity,queue,avg_waiting_time_in_min,congestion_duration_in_h,severe_congestion_duration_in_h,VMT,VHT,PMT,PHT,PDT_vf,PDT_vc,geometry,");

//        travel_time_per_iteration_map[iteration_k]

        //ODME
        if (assignment.assignment_mode == 3)
            fprintf(g_pFileLinkMOE, "obs_count,upper_type,dev,");

        for (int iteration_number = 0; iteration_number < min(20,assignment.g_number_of_column_generation_iterations); iteration_number++)
            fprintf(g_pFileLinkMOE, "TT_%d,", iteration_number);

        fprintf(g_pFileLinkMOE, "notes\n");

        //Initialization for all nodes
        for (int i = 0; i < g_link_vector.size(); ++i)
        {
            // virtual connectors
            if (g_link_vector[i].link_type == -1)
                continue;

            for (int tau = 0; tau < assignment.g_number_of_demand_periods; ++tau)
            {
                float speed = g_link_vector[i].free_speed;  // default speed 

                if (g_link_vector[i].VDF_period[tau].avg_travel_time > 0.001f)
                    speed = g_link_vector[i].link_distance_VDF / (g_link_vector[i].VDF_period[tau].avg_travel_time / 60.0);

                float speed_ratio = speed / max(1, g_link_vector[i].free_speed);  // default speed 
                float vehicle_volume = g_link_vector[i].vehicle_flow_volume_per_period[tau] + g_link_vector[i].VDF_period[tau].preload;
                float person_volume = g_link_vector[i].person_flow_volume_per_period[tau] + g_link_vector[i].VDF_period[tau].preload;
                //VMT,VHT,PMT,PHT,PDT
                float VMT = vehicle_volume * g_link_vector[i].link_distance_VDF;

                float VHT = vehicle_volume * g_link_vector[i].VDF_period[tau].avg_travel_time / 60.0;
                float PMT = person_volume * g_link_vector[i].link_distance_VDF;
                float PHT = person_volume * g_link_vector[i].VDF_period[tau].avg_travel_time / 60.0;
                float PDT_vf = person_volume * (g_link_vector[i].VDF_period[tau].avg_travel_time- g_link_vector[i].VDF_period[tau].FFTT) / 60.0;
                float PDT_vc = person_volume * (g_link_vector[i].VDF_period[tau].avg_travel_time - g_link_vector[i].VDF_period[tau].VCTT) / 60.0;

                fprintf(g_pFileLinkMOE, "%s,%d,%d,%s,%s,%s,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.1f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,\"%s\",",
                    g_link_vector[i].link_id.c_str(),
                    g_node_vector[g_link_vector[i].from_node_seq_no].node_id,
                    g_node_vector[g_link_vector[i].to_node_seq_no].node_id,
                    g_link_vector[i].link_type_name.c_str(),
                    g_link_vector[i].link_type_code.c_str(),
                    assignment.g_DemandPeriodVector[tau].time_period.c_str(),
                    vehicle_volume,
                    g_link_vector[i].person_flow_volume_per_period[tau] + g_link_vector[i].VDF_period[tau].preload,
                    g_link_vector[i].VDF_period[tau].avg_travel_time,
                    speed ,  /* 60.0 is used to convert min to hour */
                    speed_ratio,
                    g_link_vector[i].VDF_period[tau].VOC,
                    g_link_vector[i].VDF_period[tau].DOC,
                    g_link_vector[i].VDF_period[tau].lane_based_ultimate_hourly_capacity,
                    g_link_vector[i].VDF_period[tau].queue_length,
                    g_link_vector[i].VDF_period[tau].avg_waiting_time,
                    g_link_vector[i].VDF_period[tau].P,
                    g_link_vector[i].VDF_period[tau].Severe_Congestion_P,

                    VMT, VHT, PMT, PHT, PDT_vf, PDT_vc,
                    g_link_vector[i].geometry.c_str());

                if (assignment.assignment_mode == 3)  //ODME
                {
                    if (g_link_vector[i].obs_count >= 1) //ODME
                        fprintf(g_pFileLinkMOE, "%.1f,%d,%.1f,", g_link_vector[i].obs_count, g_link_vector[i].upper_bound_flag, g_link_vector[i].est_count_dev);
                    else
                        fprintf(g_pFileLinkMOE, ",,");
                }

                for (int iteration_number = 0; iteration_number < min(20, assignment.g_number_of_column_generation_iterations); iteration_number++)
                {
                    float tt = -1;
                    if(g_link_vector[i].VDF_period[tau].travel_time_per_iteration_map.find(iteration_number) != g_link_vector[i].VDF_period[tau].travel_time_per_iteration_map.end())
                    { 
                        tt = g_link_vector[i].VDF_period[tau].travel_time_per_iteration_map[iteration_number];
                    }
                    else
                    {
                        tt = -1;
                    }

                    fprintf(g_pFileLinkMOE, "%.2f,", tt);
                }
                    

                fprintf(g_pFileLinkMOE, "period-based\n");
            }

        }
        fclose(g_pFileLinkMOE);
    }

    if (assignment.assignment_mode == 0 || assignment.path_output == 0)  //LUE
    {
        FILE* g_pFilePathMOE = nullptr;
        fopen_ss(&g_pFilePathMOE, "route_assignment.csv", "w");
        fopen_ss(&g_pFilePathMOE, "no data under lue mode", "w");
        fclose(g_pFilePathMOE);

    }
    else if (assignment.assignment_mode >= 1)  //UE mode, or ODME, DTA
    {
        dtalog.output() << "writing route_assignment.csv.." << endl;

        double path_time_vector[MAX_LINK_SIZE_IN_A_PATH];
        float path_node_vector[MAX_LINK_SIZE_IN_A_PATH];
        FILE* g_pFilePathMOE = nullptr;
        fopen_ss(&g_pFilePathMOE, "route_assignment.csv", "w");

        if (!g_pFilePathMOE)
        {
            dtalog.output() << "File route_assignment.csv cannot be opened." << endl;
            g_program_stop();
        }

        fprintf(g_pFilePathMOE, "first_column,path_no,o_zone_id,d_zone_id,within_OD_path_no,path_id,information_type,agent_type,demand_period,volume,assign_volume,subarea_flag,sensor_flag,toll,#_of_nodes,travel_time,VDF_travel_time,VDF_travel_time_without_access_link,distance_km,node_sequence,link_sequence,");

        for (int iteration_number = 0; iteration_number < min(20, assignment.g_number_of_column_updating_iterations); iteration_number++)
        { 
            fprintf(g_pFilePathMOE, "TT_%d,", iteration_number);
            fprintf(g_pFilePathMOE, "Vol_%d,", iteration_number);
        }

        fprintf(g_pFilePathMOE, "geometry,");

        fprintf(g_pFilePathMOE, "link_type_name_sequence,link_code_sequence,link_link_distance_VDF_sequence,link_FFTT_sequence,");
        fprintf(g_pFilePathMOE, "\n");

        int count = 1;

        clock_t start_t, end_t;
        start_t = clock();
        clock_t iteration_t;


        int agent_type_size = assignment.g_AgentTypeVector.size();
        int zone_size = g_zone_vector.size();
        int demand_period_size = assignment.g_DemandPeriodVector.size();

        CColumnVector* p_column_pool;

        float path_toll = 0;
        float path_distance = 0;
        float path_travel_time = 0;
        float path_travel_time_without_access_link = 0;
        float path_delay = 0;
        float path_FF_travel_time = 0;
        float time_stamp = 0;

        std::map<int, CColumnPath>::iterator it, it_begin, it_end;

        if (assignment.major_path_volume_threshold > 0.00001)  // performing screening of path flow pattern
        {

            //initialization 
            bool b_subarea_mode = false;

            int number_of_links = g_link_vector.size();
            for (int i = 0; i < number_of_links; ++i)
            {
                for (int tau = 0; tau < demand_period_size; ++tau)
                {
                    // used in travel time calculation
                    g_link_vector[i].background_vehicle_flow_volume_per_period[tau] = 0;
                }

                if (g_node_vector[g_link_vector[i].from_node_seq_no].subarea_id >= 1 && g_node_vector[g_link_vector[i].to_node_seq_no].node_id >= 1)
                {

                    g_link_vector[i].subarea_id = g_node_vector[g_link_vector[i].from_node_seq_no].subarea_id;
                    b_subarea_mode = true;
                }
                else
                    g_link_vector[i].subarea_id = 0;

            }


            /// <summary>  screening the path flow pattern
            for (int orig = 0; orig < zone_size; ++orig)
            {

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
                                    int subarea_output_flag = 0;
                                    if (b_subarea_mode == true)
                                    {
                                        int insubarea_flag = 0;

                                        for (int nl = 0; nl < it->second.m_link_size; ++nl)  // arc a
                                        {
                                            int link_seq_no = it->second.path_link_vector[nl];

                                            if (g_link_vector[link_seq_no].subarea_id >= 1)
                                            {
                                                insubarea_flag = 1;
                                                break;
                                            }

                                        }
                                        // 
                                        if (insubarea_flag && it->second.path_volume > assignment.major_path_volume_threshold)
                                        {
                                            subarea_output_flag = 1;
                                        }

                                    }
                                    else
                                    {
                                        if (it->second.path_volume > assignment.major_path_volume_threshold)
                                            subarea_output_flag = 1;

                                    }
                                    if (subarea_output_flag == 0)
                                    {
                                        it->second.subarea_output_flag = 0;  // disable the output of this column into route_assignment.csv

                                        for (int nl = 0; nl < it->second.m_link_size; ++nl)  // arc a
                                        {
                                            int link_seq_no = it->second.path_link_vector[nl];
                                            g_link_vector[link_seq_no].background_vehicle_flow_volume_per_period[tau] += it->second.path_volume;
                                        }
                                    }

                                }
                            }
                        }
                        /// </summary>
                        /// <param name="assignment"></param>
                    }

                }
            }

            /// output background_link_volume.csv
            dtalog.output() << "writing link_performance.csv.." << endl;

            int b_debug_detail_flag = 0;
            FILE* g_pFileLinkMOE = nullptr;
            int b_background_link_volume_file = 0;

            if(b_background_link_volume_file)
            {
            fopen_ss(&g_pFileLinkMOE, "link_background_volume.csv", "w");
            if (!g_pFileLinkMOE)
            {
                dtalog.output() << "File link_background_volume.csv cannot be opened." << endl;
                g_program_stop();
            }
            else
            {
                fprintf(g_pFileLinkMOE, "link_id,from_node_id,to_node_id,from_cell_code,time_period,volume,background_volume,major_path_volume,ratio_of_major_path_flow,geometry,");

                fprintf(g_pFileLinkMOE, "notes\n");

                //Initialization for all nodes
                for (int i = 0; i < g_link_vector.size(); ++i)
                {
                  //   virtual connectors
                    if (g_link_vector[i].link_type == -1)
                        continue;

                    for (int tau = 0; tau < assignment.g_number_of_demand_periods; ++tau)
                    {
                        double volume = g_link_vector[i].vehicle_flow_volume_per_period[tau] + g_link_vector[i].VDF_period[tau].preload;
                        double major_path_link_volume = g_link_vector[i].vehicle_flow_volume_per_period[tau] + g_link_vector[i].VDF_period[tau].preload - g_link_vector[i].background_vehicle_flow_volume_per_period[tau];
                        double ratio = major_path_link_volume / max(volume, 0.000001);

                        if (volume < 0.0000001)
                            ratio = -1;
                        fprintf(g_pFileLinkMOE, "%s,%d,%d,%s,%s,%.3f,%.3f,%.3f,%.3f,\"%s\",",
                            g_link_vector[i].link_id.c_str(),
                            g_node_vector[g_link_vector[i].from_node_seq_no].node_id,
                            g_node_vector[g_link_vector[i].to_node_seq_no].node_id,
                            g_node_vector[g_link_vector[i].from_node_seq_no].cell_str.c_str(),
                            assignment.g_DemandPeriodVector[tau].time_period.c_str(),
                            g_link_vector[i].vehicle_flow_volume_per_period[tau] + g_link_vector[i].VDF_period[tau].preload,
                            g_link_vector[i].background_vehicle_flow_volume_per_period[tau],
                            major_path_link_volume,
                            ratio,
                            g_link_vector[i].geometry.c_str());
                        fprintf(g_pFileLinkMOE, "\n");

                    }

                }

                fclose(g_pFileLinkMOE);
            }

            }


        } // end of path flow pattern screening 
        dtalog.output() << "writing data for " << zone_size << "  zones " << endl;

        for (int orig = 0; orig < zone_size; ++orig)
        {
            if (g_zone_vector[orig].zone_id % 100 == 0)
                dtalog.output() << "o zone id =  " << g_zone_vector[orig].zone_id << endl;

            for (int at = 0; at < agent_type_size; ++at)
            {
                for (int dest = 0; dest < zone_size; ++dest)
                {
                    for (int tau = 0; tau < demand_period_size; ++tau)
                    {
                        p_column_pool = &(assignment.g_column_pool[orig][dest][at][tau]);
                        if (p_column_pool->od_volume > 0 ||
                            assignment.zone_seq_no_2_activity_mapping.find(orig) != assignment.zone_seq_no_2_activity_mapping.end() ||
                            assignment.zone_seq_no_2_info_mapping.find(dest) != assignment.zone_seq_no_2_info_mapping.end() || 
                                (assignment.zone_seq_no_2_info_mapping.find(orig) != assignment.zone_seq_no_2_info_mapping.end()
                                    && assignment.g_AgentTypeVector[at].real_time_information >= 1)
                            )
                                
                        {

                            int information_type = p_column_pool->information_type;


                            time_stamp = (assignment.g_DemandPeriodVector[tau].starting_time_slot_no + assignment.g_DemandPeriodVector[tau].ending_time_slot_no) / 2.0 * MIN_PER_TIMESLOT;

                            // scan through the map with different node sum for different continuous paths
                            it_begin = p_column_pool->path_node_sequence_map.begin();
                            it_end = p_column_pool->path_node_sequence_map.end();

                            for (it = it_begin; it != it_end; ++it)
                            {
                                if (it->second.subarea_output_flag == 0)
                                    continue;

                                if (count % 100000 == 0)
                                {
                                    end_t = clock();
                                    iteration_t = end_t - start_t;
                                    dtalog.output() << "writing " << count / 1000 << "K agents with CPU time " << iteration_t / 1000.0 << " s" << endl;
                                }

                                path_toll = 0;
                                path_distance = 0;
                                path_travel_time = 0;
                                path_travel_time_without_access_link = 0;
                                path_FF_travel_time = 0;

                                path_time_vector[0] = time_stamp;
                                path_time_vector[1] = time_stamp;

                                string link_code_str;


                                if (it->second.m_link_size >= MAX_LINK_SIZE_IN_A_PATH-2)
                                {
                                    dtalog.output() << "error: it->second.m_link_size < MAX_LINK_SIZE_IN_A_PATH" << endl;
                                    dtalog.output() << "o= " << g_zone_vector[orig].zone_id << endl;
                                    dtalog.output() << "d= " << g_zone_vector[dest].zone_id << endl;
                                    dtalog.output() << "agent type= " << assignment.g_AgentTypeVector[at].agent_type.c_str() << endl;

                                    for (int nl = 0; nl < it->second.m_link_size; ++nl)  // arc a
                                    {
                                        dtalog.output() << "node no." << nl << " =" << g_node_vector[it->second.path_node_vector[nl]].node_id << endl;
                                    }

                                    g_program_stop();
                                }

                                for (int nl = 0; nl < it->second.m_link_size; ++nl)  // arc a
                                {

                                    int link_seq_no = it->second.path_link_vector[nl];
                                    if (g_link_vector[link_seq_no].link_type >= 0)
                                    {
                                        path_toll += g_link_vector[link_seq_no].VDF_period[tau].toll[at];
                                        path_distance += g_link_vector[link_seq_no].link_distance_VDF;
                                        float link_travel_time = g_link_vector[link_seq_no].travel_time_per_period[tau];
                                        path_travel_time += link_travel_time;
                                        
                                        if(g_link_vector[link_seq_no].link_type !=1000)   // skip access link
                                        {
                                            path_travel_time_without_access_link += link_travel_time;
                                        }

                                        path_FF_travel_time += g_link_vector[link_seq_no].VDF_period[tau].FFTT;
                                        time_stamp += link_travel_time;
                                        path_time_vector[nl+1] = time_stamp;
                                        link_code_str += g_link_vector[link_seq_no].link_code_str;
                                    }
                                    else
                                    {
                                        int virtual_link = 0;
                                    }
                                }

                                double total_agent_path_travel_time = 0;

                                for (int vi = 0; vi < it->second.agent_simu_id_vector.size(); ++vi)
                                {
                                    int agent_simu_id = it->second.agent_simu_id_vector[vi];
                                    CAgent_Simu* pAgentSimu = g_agent_simu_vector[agent_simu_id];
                                    total_agent_path_travel_time += pAgentSimu->path_travel_time_in_min;
                                }

                                double final_path_travel_time = path_travel_time;  // by default

                                if (it->second.agent_simu_id_vector.size() > 1)  // with simulated agents
                                {
                                    final_path_travel_time = total_agent_path_travel_time / it->second.agent_simu_id_vector.size();
                                }



                                int virtual_first_link_delta = 1;
                                int virtual_last_link_delta = 1;
                                // fixed routes have physical nodes always, without virtual connectors
                                if (p_column_pool->bfixed_route)
                                {

                                    virtual_first_link_delta = 0;
                                    virtual_last_link_delta = 1;
                                }
                                if (information_type == 1 && it->second.path_volume < 0.1)
                                    continue;

                                if (information_type == 1)  // information diversion start from physical nodes
                                {
                                    virtual_first_link_delta = 0;
                                    virtual_last_link_delta = 1;

                                    it->second.path_volume = it->second.diverted_vehicle_map.size();
                                }

                                if (p_column_pool->activity_zone_no_vector.size()>0)  // with activity zones 
                                {
                                    virtual_first_link_delta = 0;
                                    virtual_last_link_delta = 0;
                                }


                                if (assignment.zone_seq_no_2_info_mapping.find(orig) != assignment.zone_seq_no_2_info_mapping.end()
                                    && assignment.g_AgentTypeVector[at].real_time_information >= 1)
                                {
                                    if(it->second.path_volume <1)
                                        it->second.path_volume = 1;  // reset the volume as 1 to enable visualization 

                                }

                                // assignment_mode = 1, path flow mode
                                if (assignment.assignment_mode >= 1)
                                {
                                    if (it->second.m_node_size - virtual_first_link_delta - virtual_last_link_delta <= 2)
                                        continue;

                                    fprintf(g_pFilePathMOE, ",%d,%d,%d,%d,%s,%d,%s,%s,%.2f,%d,%d,%d,%.1f,%d,%.1f,%.4f,%.4f,%.4f,",
                                        count,
                                        g_zone_vector[orig].zone_id,
                                        g_zone_vector[dest].zone_id,
                                        it->second.path_seq_no,
                                        it->second.path_id.c_str(),
                                        information_type,
                                        assignment.g_AgentTypeVector[at].agent_type.c_str(),
                                        assignment.g_DemandPeriodVector[tau].demand_period.c_str(),
                                        it->second.path_volume,
                                        it->second.agent_simu_id_vector.size(),
                                        it->second.subarea_output_flag,
                                        it->second.measurement_flag,
                                        path_toll,
                                        it->second.m_node_size- virtual_first_link_delta- virtual_last_link_delta,
                                        final_path_travel_time,
                                        path_travel_time,
                                        path_travel_time_without_access_link,
                                        //path_FF_travel_time,
                                        //final_path_travel_time- path_FF_travel_time,
                                        path_distance);

                                    /* Format and print various data */
                                    for (int ni = 0 + virtual_first_link_delta; ni < it->second.m_node_size - virtual_last_link_delta; ++ni)
                                        fprintf(g_pFilePathMOE, "%d;", g_node_vector[it->second.path_node_vector[ni]].node_id);

                                    fprintf(g_pFilePathMOE, ",");
                                    int link_seq_no;
                                    // link id sequence
                                    for (int nl = 0 + virtual_first_link_delta; nl < it->second.m_link_size - virtual_last_link_delta; ++nl)
                                    {
                                        link_seq_no = it->second.path_link_vector[nl];
                                        fprintf(g_pFilePathMOE, "%s;", g_link_vector[link_seq_no].link_id.c_str());
                                    }
                                    fprintf(g_pFilePathMOE, ",");


                                    //for (int nt = 0 + virtual_first_link_delta; nt < min(5000-1, it->second.m_node_size - virtual_last_link_delta); ++nt)
                                    //{
                                    //    if(path_time_vector[nt] < assignment.g_LoadingEndTimeInMin)
                                    //    {
                                    //    fprintf(g_pFilePathMOE, "%s;", g_time_coding(path_time_vector[nt]).c_str());
                                    //    }
                                    //}
                                    //fprintf(g_pFilePathMOE, ",");

                                    //for (int nt = 0 + virtual_first_link_delta; nt < it->second.m_link_size+1 - virtual_last_link_delta; ++nt)
                                    //    fprintf(g_pFilePathMOE, "%.2f;", path_time_vector[nt]);

                                    //fprintf(g_pFilePathMOE, ",");

                                    //for (int nt = 0 + virtual_first_link_delta; nt < it->second.m_link_size - virtual_last_link_delta; ++nt)
                                    //    fprintf(g_pFilePathMOE, "%.2f;", path_time_vector[nt+1]- path_time_vector[nt]);

                                    //fprintf(g_pFilePathMOE, ",");

                                    for (int iteration_number = 0; iteration_number < min(20, assignment.g_number_of_column_updating_iterations); iteration_number++)
                                    {
                                        double TT = -1;
                                        double Vol = 0;
                                        if (it->second.path_gradient_cost_per_iteration_map.find(iteration_number) != it->second.path_gradient_cost_per_iteration_map.end())
                                        {
                                            TT = it->second.path_gradient_cost_per_iteration_map[iteration_number];
                                        }

                                        if (it->second.path_volume_per_iteration_map.find(iteration_number) != it->second.path_volume_per_iteration_map.end())
                                        {
                                            Vol = it->second.path_volume_per_iteration_map[iteration_number];
                                        }

                                        fprintf(g_pFilePathMOE, "%f,%f,", TT, Vol);

                                    }

                                    if(it->second.m_node_size- virtual_first_link_delta- virtual_last_link_delta>=2)
                                    {
                                    fprintf(g_pFilePathMOE, "\"LINESTRING (");
                                        for (int ni = 0 + virtual_first_link_delta; ni < it->second.m_node_size - virtual_last_link_delta; ++ni)
                                        {
                                            fprintf(g_pFilePathMOE, "%f %f", g_node_vector[it->second.path_node_vector[ni]].x,
                                                g_node_vector[it->second.path_node_vector[ni]].y);

                                            if (ni != it->second.m_node_size - virtual_last_link_delta - 1)
                                                fprintf(g_pFilePathMOE, ", ");
                                        }

                                    fprintf(g_pFilePathMOE, ")\"");
                                    }

                                    if(it->second.path_volume>=5)  // critcal path volume
                                    {
                                        // link type name sequenece
                                        for (int nl = 0 + virtual_first_link_delta; nl < it->second.m_link_size - virtual_last_link_delta; ++nl)
                                        {
                                            link_seq_no = it->second.path_link_vector[nl];
                                            fprintf(g_pFilePathMOE, "%s;", g_link_vector[link_seq_no].link_type_name.c_str());
                                        }
                                        fprintf(g_pFilePathMOE, ",");

                                        // link code sequenece
                                        for (int nl = 0 + virtual_first_link_delta; nl < it->second.m_link_size - virtual_last_link_delta; ++nl)
                                        {
                                            link_seq_no = it->second.path_link_vector[nl];
                                            fprintf(g_pFilePathMOE, "%s;", g_link_vector[link_seq_no].link_code_str.c_str());
                                        }
                                        fprintf(g_pFilePathMOE, ",");

                                        // link link_distance_VDF sequenece
                                        for (int nl = 0 + virtual_first_link_delta; nl < it->second.m_link_size - virtual_last_link_delta; ++nl)
                                        {
                                            link_seq_no = it->second.path_link_vector[nl];
                                            fprintf(g_pFilePathMOE, "%.3f;", g_link_vector[link_seq_no].link_distance_VDF);
                                        }
                                        fprintf(g_pFilePathMOE, ",");

                                        // link FFTT sequenece
                                        for (int nl = 0 + virtual_first_link_delta; nl < it->second.m_link_size - virtual_last_link_delta; ++nl)
                                        {
                                            link_seq_no = it->second.path_link_vector[nl];
                                            fprintf(g_pFilePathMOE, "%.3f;", g_link_vector[link_seq_no].free_flow_travel_time_in_min);
                                        }
                                        fprintf(g_pFilePathMOE, ",");
                                    }
                                    fprintf(g_pFilePathMOE, "\n");
                                    count++;
                                }

                            }
                        }
                    }
                }
            }
        }
        fclose(g_pFilePathMOE);
    }

    if (assignment.VDF_type == 2)  //numerical queue evolution mode
        g_output_dynamic_queue_profile();
}

void g_output_accessibility_result(Assignment& assignment)
{

    if (assignment.assignment_mode == 0 || assignment.path_output == 0)  //LUE
    {
        FILE* g_pFilePathMOE = nullptr;
        fopen_ss(&g_pFilePathMOE, "od_accessibility.csv", "w");
        fopen_ss(&g_pFilePathMOE, "no data under lue mode", "w");
        fclose(g_pFilePathMOE);

    }
    else if (assignment.assignment_mode >= 1)  //UE mode, or ODME, DTA
    {
        dtalog.output() << "writing od_accessibility.csv.." << endl;

        double path_time_vector[MAX_LINK_SIZE_IN_A_PATH];
        FILE* g_pFilePathMOE = nullptr;
        fopen_ss(&g_pFilePathMOE, "od_accessibility.csv", "w");

        if (!g_pFilePathMOE)
        {
            dtalog.output() << "File od_accessibility.csv cannot be opened." << endl;
            g_program_stop();
        }

        fprintf(g_pFilePathMOE, "od_no,o_zone_id,d_zone_id,o_zone_agent_type_cover_flag,d_zone_agent_type_cover_flag,path_id,information_type,agent_type,demand_period,volume,assign_volume,subarea_flag,sensor_flag,toll,travel_time,VDF_travel_time,travel_time_without_access_link,distance_km\n");

        int count = 1;

        clock_t start_t, end_t;
        start_t = clock();
        clock_t iteration_t;


        int agent_type_size = assignment.g_AgentTypeVector.size();
        int zone_size = g_zone_vector.size();
        int demand_period_size = assignment.g_DemandPeriodVector.size();

        CColumnVector* p_column_pool;

        float path_toll = 0;
        float path_distance = 0;
        float path_travel_time = 0;
        float path_travel_time_without_access_link = 0;
        float path_delay = 0;
        float path_FF_travel_time = 0;
        float time_stamp = 0;

        std::map<int, CColumnPath>::iterator it, it_begin, it_end;

        if (assignment.major_path_volume_threshold > 0.00001)  // performing screening of path flow pattern
        {

            //initialization 
            bool b_subarea_mode = false;

            int number_of_links = g_link_vector.size();
            for (int i = 0; i < number_of_links; ++i)
            {
                for (int tau = 0; tau < demand_period_size; ++tau)
                {
                    // used in travel time calculation
                    g_link_vector[i].background_vehicle_flow_volume_per_period[tau] = 0;
                }

                if (g_node_vector[g_link_vector[i].from_node_seq_no].subarea_id >= 1 && g_node_vector[g_link_vector[i].to_node_seq_no].node_id >= 1)
                {

                    g_link_vector[i].subarea_id = g_node_vector[g_link_vector[i].from_node_seq_no].subarea_id;
                    b_subarea_mode = true;
                }
                else
                    g_link_vector[i].subarea_id = 0;

            }


            /// <summary>  screening the path flow pattern
            for (int orig = 0; orig < zone_size; ++orig)
            {

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
                                    int subarea_output_flag = 0;
                                    if (b_subarea_mode == true)
                                    {
                                        int insubarea_flag = 0;

                                        for (int nl = 0; nl < it->second.m_link_size; ++nl)  // arc a
                                        {
                                            int link_seq_no = it->second.path_link_vector[nl];

                                            if (g_link_vector[link_seq_no].subarea_id >= 1)
                                            {
                                                insubarea_flag = 1;
                                                break;
                                            }

                                        }
                                        // 
                                        if (insubarea_flag && it->second.path_volume > assignment.major_path_volume_threshold)
                                        {
                                            subarea_output_flag = 1;
                                        }

                                    }
                                    else
                                    {
                                        if (it->second.path_volume > assignment.major_path_volume_threshold)
                                            subarea_output_flag = 1;

                                    }
                                    if (subarea_output_flag == 0)
                                    {
                                        it->second.subarea_output_flag = 0;  // disable the output of this column into route_assignment.csv

                                        for (int nl = 0; nl < it->second.m_link_size; ++nl)  // arc a
                                        {
                                            int link_seq_no = it->second.path_link_vector[nl];
                                            g_link_vector[link_seq_no].background_vehicle_flow_volume_per_period[tau] += it->second.path_volume;
                                        }
                                    }

                                }
                            }
                        }
                        /// </summary>
                        /// <param name="assignment"></param>
                    }

                }
            }


        } // end of path flow pattern screening 
        dtalog.output() << "writing data for " << zone_size << "  zones " << endl;

        for (int orig = 0; orig < zone_size; ++orig)
        {
            if (g_zone_vector[orig].zone_id % 100 == 0)
                dtalog.output() << "o zone id =  " << g_zone_vector[orig].zone_id << endl;

            for (int at = 0; at < agent_type_size; ++at)
            {
                for (int dest = 0; dest < zone_size; ++dest)
                {
                    for (int tau = 0; tau < demand_period_size; ++tau)
                    {
                        p_column_pool = &(assignment.g_column_pool[orig][dest][at][tau]);
                        if (p_column_pool->od_volume > 0 || assignment.zone_seq_no_2_info_mapping.find(orig) != assignment.zone_seq_no_2_info_mapping.end())
                        {

                            int information_type = 0;

                            if (assignment.zone_seq_no_2_info_mapping.find(orig) != assignment.zone_seq_no_2_info_mapping.end())
                            {
                                information_type = 1;
                            }

                            time_stamp = (assignment.g_DemandPeriodVector[tau].starting_time_slot_no + assignment.g_DemandPeriodVector[tau].ending_time_slot_no) / 2.0 * MIN_PER_TIMESLOT;

                            // scan through the map with different node sum for different continuous paths
                            it_begin = p_column_pool->path_node_sequence_map.begin();
                            it_end = p_column_pool->path_node_sequence_map.end();

                            for (it = it_begin; it != it_end; ++it)
                            {
                                if (it->second.subarea_output_flag == 0)
                                    continue;

                                if (count % 100000 == 0)
                                {
                                    end_t = clock();
                                    iteration_t = end_t - start_t;
                                    dtalog.output() << "writing " << count / 1000 << "K agents with CPU time " << iteration_t / 1000.0 << " s" << endl;
                                }

                                path_toll = 0;
                                path_distance = 0;
                                path_travel_time = 0;
                                path_travel_time_without_access_link = 0;
                                path_FF_travel_time = 0;

                                path_time_vector[0] = time_stamp;
                                string link_code_str;

                                for (int nl = 0; nl < it->second.m_link_size; ++nl)  // arc a
                                {

                                    int link_seq_no = it->second.path_link_vector[nl];
                                    if (g_link_vector[link_seq_no].link_type >= 0)
                                    {
                                        path_toll += g_link_vector[link_seq_no].VDF_period[tau].toll[at];
                                        path_distance += g_link_vector[link_seq_no].link_distance_VDF;
                                        float link_travel_time = g_link_vector[link_seq_no].travel_time_per_period[tau];
                                        path_travel_time += link_travel_time;

                                        if (g_link_vector[link_seq_no].link_type != 1000)   // skip access link
                                        {
                                            path_travel_time_without_access_link += link_travel_time;
                                        }
                                        path_FF_travel_time += g_link_vector[link_seq_no].VDF_period[tau].FFTT;
                                        time_stamp += link_travel_time;
                                        path_time_vector[nl + 1] = time_stamp;

                                        link_code_str += g_link_vector[link_seq_no].link_code_str;
                                    }
                                    else
                                    {
                                        int virtual_link = 0;
                                    }
                                }

                                double total_agent_path_travel_time = 0;

                                for (int vi = 0; vi < it->second.agent_simu_id_vector.size(); ++vi)
                                {
                                    int agent_simu_id = it->second.agent_simu_id_vector[vi];
                                    CAgent_Simu* pAgentSimu = g_agent_simu_vector[agent_simu_id];
                                    total_agent_path_travel_time += pAgentSimu->path_travel_time_in_min;
                                }

                                double final_path_travel_time = path_travel_time;  // by default

                                if (it->second.agent_simu_id_vector.size() > 1)  // with simulated agents
                                {
                                    final_path_travel_time = total_agent_path_travel_time / it->second.agent_simu_id_vector.size();
                                }



                                int virtual_first_link_delta = 1;
                                int virtual_last_link_delta = 1;
                                // fixed routes have physical nodes always, without virtual connectors
                                if (p_column_pool->bfixed_route)
                                {

                                    virtual_first_link_delta = 0;
                                    virtual_last_link_delta = 1;
                                }
                                if (information_type == 1 && it->second.path_volume < 0.1)
                                    continue;

                                if (information_type == 1)  // information diversion start from physical nodes
                                {
                                    virtual_first_link_delta = 0;
                                    virtual_last_link_delta = 1;

                                    it->second.path_volume = it->second.diverted_vehicle_map.size();
                                }


                                // assignment_mode = 1, path flow mode
                                if (assignment.assignment_mode >= 1)
                                {

                                    int o_zone_agent_type_cover_flag = 0;
                                    int d_zone_agent_type_cover_flag = 0;

                                    if (assignment.g_AgentTypeVector[at].zone_id_cover_map.find(g_zone_vector[orig].zone_id) != assignment.g_AgentTypeVector[at].zone_id_cover_map.end())
                                    {
                                        o_zone_agent_type_cover_flag = 1;
                                    }

                                    if (assignment.g_AgentTypeVector[at].zone_id_cover_map.find(g_zone_vector[dest].zone_id) != assignment.g_AgentTypeVector[at].zone_id_cover_map.end())
                                    {
                                        d_zone_agent_type_cover_flag = 1;
                                    }


                                    fprintf(g_pFilePathMOE, "%d,%d,%d,%d,%d,%d,%d,%s,%s,%.2f,%d,%d,%d,%.1f,%.1f,%.4f,%.4f,%.4f\n",
                                        count,
                                        g_zone_vector[orig].zone_id,
                                        g_zone_vector[dest].zone_id,
                                        o_zone_agent_type_cover_flag,
                                        d_zone_agent_type_cover_flag,
                                        it->second.path_seq_no,
                                        information_type,
                                        assignment.g_AgentTypeVector[at].agent_type.c_str(),
                                        assignment.g_DemandPeriodVector[tau].demand_period.c_str(),
                                        it->second.path_volume,
                                        it->second.agent_simu_id_vector.size(),
                                        it->second.subarea_output_flag,
                                        it->second.measurement_flag,
                                        path_toll,
                                        final_path_travel_time,
                                        path_travel_time,
                                        path_travel_time_without_access_link,
                                        //path_FF_travel_time,
                                        //final_path_travel_time- path_FF_travel_time,
                                        path_distance);


                                    count++;
                                }

                            }
                        }
                    }
                }
            }
        }
        fclose(g_pFilePathMOE);
    }
}

void g_output_dynamic_link_performance(Assignment& assignment, int output_mode = 1)
{
    dtalog.output() << "writing dynamic_link_performance.csv.." << endl;

    int b_debug_detail_flag = 0;
    FILE* g_pFileLinkMOE = nullptr;

    string file_name = "dynamic_link_performance.csv";

    if (output_mode == 0)
    {
        file_name = "dynamic_link_performance_hd.csv";
    }

    if (output_mode == 2)
    {
        file_name = "dynamic_link_performance_horizon.csv";
    }

 
    fopen_ss(&g_pFileLinkMOE, file_name.c_str(), "w");

    if (!g_pFileLinkMOE)
    {
        dtalog.output() << "File " << file_name.c_str() << " cannot be opened." << endl;
        g_program_stop();
    }
    else
    {

        // Option 2: BPR-X function
        fprintf(g_pFileLinkMOE, "link_id,tmc_corridor_name,link_type_name,from_node_id,to_node_id,from_cell_code,lanes,link_distance_VDF,free_speed,FFTT,time_period,volume,CA,CD,density,queue_link_distance_VDF,discharge_rate_per_period,travel_time,waiting_time_in_min,RT_speed,speed,geometry,");
        fprintf(g_pFileLinkMOE, "notes\n");


        //Initialization for all nodes
        for (int i = 0; i < g_link_vector.size(); ++i)
        {
            // virtual connectors
            if (g_link_vector[i].link_type == -1)
                continue;

            // first loop for time t
            for (int t = 1; t < assignment.g_number_of_intervals_in_min; ++t)
            {
                if (assignment.dynamic_link_performance_sampling_interval_in_min < 1)
                    assignment.dynamic_link_performance_sampling_interval_in_min = 1;

                int sampling_time_interval = assignment.dynamic_link_performance_sampling_interval_in_min;

                if (output_mode == 0)
                    sampling_time_interval = assignment.dynamic_link_performance_sampling_interval_hd_in_min * 60; // min by min

                if (output_mode == 2)
                    sampling_time_interval = assignment.g_number_of_loading_intervals_in_sec / 60; // simulation horizon

                if (t % (sampling_time_interval) == 0)
                {
                    int time_in_min = t;  //relative time

                    float volume = 0;
                    float queue = 0;
                    float waiting_time_in_min = 0;
                    int arrival_flow_volume = 0;
                    float avg_waiting_time_in_min = 0;

                    float travel_time = (float)(g_link_vector[i].free_flow_travel_time_in_min + avg_waiting_time_in_min);
                    float speed = g_link_vector[i].link_distance_VDF / (g_link_vector[i].free_flow_travel_time_in_min / 60.0);
                    float virtual_arrival = 0;

                    float discharge_rate = g_link_vector[i].lane_capacity * g_link_vector[i].number_of_lanes * assignment.dynamic_link_performance_sampling_interval_in_min / 60.0;

                    if (time_in_min >= 1)
                    {
                        volume = assignment.m_LinkCumulativeDepartureVector[i][t] - assignment.m_LinkCumulativeDepartureVector[i][t - sampling_time_interval];


                        queue = assignment.m_LinkCumulativeArrivalVector[i][t] - assignment.m_LinkCumulativeDepartureVector[i][t];
                        //							waiting_time_in_min = queue / (max(1, volume));

                        float waiting_time_count = 0;

                        waiting_time_in_min = assignment.m_LinkTDWaitingTime[i][t] * number_of_seconds_per_interval;

                        if (output_mode == 2)
                        {
                            waiting_time_in_min = assignment.m_LinkTotalWaitingTimeVector[i];
                            avg_waiting_time_in_min = waiting_time_in_min / max(1, arrival_flow_volume);
                        }

                        arrival_flow_volume = assignment.m_LinkCumulativeArrivalVector[i][t] - assignment.m_LinkCumulativeArrivalVector[i][t - sampling_time_interval];
                        avg_waiting_time_in_min = waiting_time_in_min / max(1, arrival_flow_volume);

                        travel_time = (float)(g_link_vector[i].free_flow_travel_time_in_min + avg_waiting_time_in_min);
                        speed = g_link_vector[i].link_distance_VDF / (max(0.00001f, travel_time / 60.0));
                    }


                    float density = (assignment.m_LinkCumulativeArrivalVector[i][t] - assignment.m_LinkCumulativeDepartureVector[i][t]) / (g_link_vector[i].link_distance_VDF * g_link_vector[i].number_of_lanes);

                    if (density >= 1000 )
                    {
                        int i_debug = 1;
                    }
                    if (output_mode == 0)
                    {
                        if (assignment.m_LinkCumulativeArrivalVector[i][t] < 1000)
                            continue; //skip
                    }

                    fprintf(g_pFileLinkMOE, "%s,%s,%s,%d,%d,%s,%d,%.3f,%.3f,%.3f,%s_%s,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,\"%s\",",
                        g_link_vector[i].link_id.c_str(),
                        g_link_vector[i].tmc_corridor_name.c_str(),
                        g_link_vector[i].link_type_name.c_str(),

                        g_node_vector[g_link_vector[i].from_node_seq_no].node_id,
                        g_node_vector[g_link_vector[i].to_node_seq_no].node_id,
                        g_node_vector[g_link_vector[i].from_node_seq_no].cell_str.c_str(),
                        g_link_vector[i].number_of_lanes,
                        g_link_vector[i].link_distance_VDF,
                        g_link_vector[i].free_speed,
                        g_link_vector[i].free_flow_travel_time_in_min,

                        g_time_coding(assignment.g_LoadingStartTimeInMin + time_in_min - sampling_time_interval).c_str(),
                        g_time_coding(assignment.g_LoadingStartTimeInMin + time_in_min).c_str(),
                        volume,
                        assignment.m_LinkCumulativeArrivalVector[i][t],
                        assignment.m_LinkCumulativeDepartureVector[i][t],
                        density,
                        queue,
                        discharge_rate,
                        travel_time,
                        avg_waiting_time_in_min,
                        speed,
                        g_link_vector[i].geometry.c_str());

                    fprintf(g_pFileLinkMOE, "simulation-based\n");
                }
            }  // for each time t
        }  // for each link l
        fclose(g_pFileLinkMOE);
    }//assignment mode 2 as simulation

}


void g_output_dynamic_link_state(Assignment& assignment, int output_mode = 1)
{
    dtalog.output() << "writing dynamic_link_state.csv.." << endl;

    int b_debug_detail_flag = 0;
    FILE* g_pFileLinkMOE = nullptr;

    string file_name = "dynamic_link_state.csv";

    fopen_ss(&g_pFileLinkMOE, file_name.c_str(), "w");

    if (!g_pFileLinkMOE)
    {
        dtalog.output() << "File " << file_name.c_str() << " cannot be opened." << endl;
        g_program_stop();
    }
    else
    {

        // Option 2: BPR-X function
        fprintf(g_pFileLinkMOE, "link_id,from_node_id,to_node_id,time_period,duration_in_sec,state,state_code\n");


        //Initialization for all nodes
        for (unsigned li = 0; li < g_link_vector.size(); ++li)
        {

            if (g_link_vector[li].timing_arc_flag || g_link_vector[li].m_link_pedefined_capacity_map.size() > 0)  // only output the capaicty for signalized data 
            {
                // reset for signalized links (not freeway links as type code != 'f' for the case of freeway workzones)
                // only for the loading period

                int t = 0;
                int last_t = t;
                int current_state = assignment.m_LinkOutFlowState[li][t];
                if (current_state == 0)  //close in the beginning
                    current_state = -1; // change the initial state to be able to be record the change below

                while (t < assignment.g_number_of_loading_intervals_in_sec - 1)
                {
                    int next_state = assignment.m_LinkOutFlowState[li][t + 1];

                    bool print_out_flag = false;

                    if (current_state == 0 && t%60==0 && g_link_vector[li].m_link_pedefined_capacity_map.find(t) != g_link_vector[li].m_link_pedefined_capacity_map.end())
                    {
                        print_out_flag = true;
                    }

                    if (print_out_flag==false && next_state == current_state && t < assignment.g_number_of_loading_intervals_in_sec - 2)
                    {
                        // do nothing 
                    }
                    else
                    {  // change of state

                        string state_str;
                        if(g_link_vector[li].timing_arc_flag)
                        { 


                        if (current_state == 1)
                            state_str = "g";

                        if (current_state == 0)
                            state_str = "r";
                        }
                        else if (g_link_vector[li].m_link_pedefined_capacity_map.find(t)!= g_link_vector[li].m_link_pedefined_capacity_map.end())
                        {
                            state_str = "w";
                        }

                        if(state_str.size()>0)  // with data to output
                        {
                        fprintf(g_pFileLinkMOE, "%s,%d,%d,%s_%s,%d,%d,%s\n",
                            g_link_vector[li].link_id.c_str(),
                            g_node_vector[g_link_vector[li].from_node_seq_no].node_id,
                            g_node_vector[g_link_vector[li].to_node_seq_no].node_id,
                            g_time_coding(assignment.g_LoadingStartTimeInMin + last_t / 60.0).c_str(),
                            g_time_coding(assignment.g_LoadingStartTimeInMin + (t + 1) / 60.0).c_str(),
                            t + 1 - last_t,
                            current_state,
                            state_str.c_str());
                        }

                        last_t = t + 1;
                        current_state = assignment.m_LinkOutFlowState[li][t + 1];

                        if (t + 1 == assignment.g_number_of_loading_intervals_in_sec - 2)
                        {
                            //boundary condition anyway
                            current_state = -1;
                        }

                    }
                    t++;
                }


            } 



        }
        fclose(g_pFileLinkMOE);

    }
}

void g_output_simulation_agents(Assignment& assignment)
{
    if (assignment.assignment_mode == 0 || assignment.trajectory_output == 0)  //LUE
    {
        FILE* g_pFilePathMOE = nullptr;
        fopen_ss(&g_pFilePathMOE, "agent.csv", "w");
        fclose(g_pFilePathMOE);
    }
    else if (assignment.assignment_mode >= 1)  //UE mode, or ODME, DTA
    {
        dtalog.output() << "writing agent.csv.." << endl;

        double path_time_vector[MAX_LINK_SIZE_IN_A_PATH];
        FILE* g_pFileAgent = nullptr;
        fopen_ss(&g_pFileAgent, "agent.csv", "w");

        if (!g_pFileAgent)
        {
            dtalog.output() << "File agent.csv cannot be opened." << endl;
            g_program_stop();
        }

        fprintf(g_pFileAgent, "agent_id,o_zone_id,d_zone_id,OD_index,path_id,#_of_links,diversion_flag,agent_type,demand_period,volume,toll,travel_time,distance_km,speed,departure_time_in_min,arrival_time_in_min,departure_time_slot_no,\n");

        int count = 1;

        clock_t start_t, end_t;
        start_t = clock();
        clock_t iteration_t;

        int agent_type_size = assignment.g_AgentTypeVector.size();
        int zone_size = g_zone_vector.size();
        int demand_period_size = assignment.g_DemandPeriodVector.size();

        CColumnVector* p_column_pool;

        float path_toll = 0;
        float path_distance = 0;
        float path_travel_time = 0;
        float path_travel_time_without_access_link = 0;
        float time_stamp = 0;

        if (assignment.trajectory_sampling_rate < 0.01)
            assignment.trajectory_sampling_rate = 0.01;
        int sampling_step = 100 / int(100 * assignment.trajectory_sampling_rate + 0.5);

        std::map<int, CColumnPath>::iterator it, it_begin, it_end;

        dtalog.output() << "writing data for " << zone_size << "  zones " << endl;

        for (int orig = 0; orig < zone_size; ++orig)
        {
            if (g_zone_vector[orig].zone_id % 100 == 0)
                dtalog.output() << "o zone id =  " << g_zone_vector[orig].zone_id << endl;

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
                                if (count % 100000 == 0)
                                {
                                    end_t = clock();
                                    iteration_t = end_t - start_t;
                                    dtalog.output() << "writing " << count / 1000 << "K agents with CPU time " << iteration_t / 1000.0 << " s" << endl;
                                }

                                if (count % sampling_step != 0)
                                    continue;


                                path_toll = 0;
                                path_distance = 0;
                                path_travel_time = 0;
                                path_travel_time_without_access_link = 0;
                                path_time_vector[0] = time_stamp;



                                // assignment_mode = 1, path flow mode
                                {
                                    // assignment_mode = 2, simulated agent flow mode //DTA simulation 

                                    for (int vi = 0; vi < it->second.agent_simu_id_vector.size(); ++vi)
                                    {
                                        int agent_simu_id = it->second.agent_simu_id_vector[vi];
                                        CAgent_Simu* pAgentSimu = g_agent_simu_vector[agent_simu_id];
                                        time_stamp = assignment.g_LoadingStartTimeInMin + pAgentSimu->departure_time_in_min;

                                        float departure_time_in_min = time_stamp;
                                        float arrival_time_in_min = assignment.g_LoadingStartTimeInMin + pAgentSimu->arrival_time_in_min;
                                        int departure_time_in_slot_no = time_stamp / MIN_PER_TIMESLOT;
                                        float speed = pAgentSimu->path_distance / max(0.001, pAgentSimu->path_travel_time_in_min) * 60;

                                        if (it->second.m_link_size >= MAX_LINK_SIZE_IN_A_PATH -1)
                                        {
                                            dtalog.output() << "error: it->second.m_link_size >= MAX_LINK_SIZE_IN_A_PATH-1" << endl;
                                            g_program_stop();
                                        }

                                        fprintf(g_pFileAgent, "%d,%d,%d,%d->%d,%d,%d,%d,%s,%s,1,%.1f,%.4f,%.4f,%.4f,%.4f,%.4f,%d",
                                            pAgentSimu->agent_id,
                                            g_zone_vector[orig].zone_id,
                                            g_zone_vector[dest].zone_id,
                                            g_zone_vector[orig].zone_id,
                                            g_zone_vector[dest].zone_id,
                                            it->second.path_seq_no,
                                            it->second.m_link_size,
                                            pAgentSimu->diversion_flag,
                                            assignment.g_AgentTypeVector[at].agent_type.c_str(),
                                            assignment.g_DemandPeriodVector[tau].demand_period.c_str(),
                                            path_toll,
                                            pAgentSimu->path_travel_time_in_min,
                                            pAgentSimu->path_distance, speed,
                                            departure_time_in_min,
                                            arrival_time_in_min,
                                            departure_time_in_slot_no);

                                        count++;

                                        fprintf(g_pFileAgent, "\n");
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        fclose(g_pFileAgent);
    }
}


void g_output_trajectory_csv(Assignment& assignment)
{

    if (assignment.assignment_mode == 0 || assignment.trajectory_output == 0)  //LUE
    {
        FILE* g_pFilePathMOE = nullptr;
        fopen_ss(&g_pFilePathMOE, "trajectory.csv", "w");
        fclose(g_pFilePathMOE);
    }
    else if (assignment.assignment_mode >= 1)  //UE mode, or ODME, DTA
    {
        dtalog.output() << "writing trajectory.csv.." << endl;

        double path_time_vector[MAX_LINK_SIZE_IN_A_PATH];
        FILE* g_pFilePathMOE = nullptr;
        fopen_ss(&g_pFilePathMOE, "trajectory.csv", "w");

        if (!g_pFilePathMOE)
        {
            dtalog.output() << "File trajectory.csv cannot be opened." << endl;
            g_program_stop();
        }

        fprintf(g_pFilePathMOE, "first_column,agent_id,o_zone_id,d_zone_id,path_id,display_code,agent_type,PCE_unit,demand_period,volume,toll,travel_time,distance_km,node_sequence,time_sequence\n");

        int count = 1;

        clock_t start_t, end_t;
        start_t = clock();
        clock_t iteration_t;

        int agent_type_size = assignment.g_AgentTypeVector.size();
        int zone_size = g_zone_vector.size();
        int demand_period_size = assignment.g_DemandPeriodVector.size();

        CColumnVector* p_column_pool;

        float path_toll = 0;
        float path_distance = 0;
        float path_travel_time = 0;
        float time_stamp = 0;

        if (assignment.trajectory_sampling_rate < 0.01)
            assignment.trajectory_sampling_rate = 0.01;
        int sampling_step = 100 / int(100 * assignment.trajectory_sampling_rate + 0.5);

        std::map<int, CColumnPath>::iterator it, it_begin, it_end;

        dtalog.output() << "writing data for " << zone_size << "  zones " << endl;

        for (int orig = 0; orig < zone_size; ++orig)
        {
            if (g_zone_vector[orig].zone_id % 100 == 0)
                dtalog.output() << "o zone id =  " << g_zone_vector[orig].zone_id << endl;

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
                                if (count % 100000 == 0)
                                {
                                    end_t = clock();
                                    iteration_t = end_t - start_t;
                                    dtalog.output() << "writing " << count / 1000 << "K agents with CPU time " << iteration_t / 1000.0 << " s" << endl;
                                }

                                if (count % sampling_step != 0)
                                    continue;


                                path_toll = 0;
                                path_distance = 0;
                                path_travel_time = 0;
                                path_time_vector[0] = time_stamp;

                                for (int nl = 0; nl < it->second.m_link_size; ++nl)  // arc a
                                {
                                    int link_seq_no = it->second.path_link_vector[nl];
                                    path_toll += g_link_vector[link_seq_no].VDF_period[tau].toll[at];
                                    path_distance += g_link_vector[link_seq_no].link_distance_VDF;
                                    float link_travel_time = g_link_vector[link_seq_no].travel_time_per_period[tau];

                                    time_stamp += link_travel_time;
                                    path_time_vector[nl + 1] = time_stamp;
                                }

                                int virtual_link_delta = 1;
                                int virtual_begin_link_delta = 1;
                                int virtual_end_link_delta = 1;
                                // fixed routes have physical nodes always, without virtual connectors
                                if (p_column_pool->bfixed_route)
                                {
                                    virtual_begin_link_delta = 0;
                                    virtual_end_link_delta = 1;

                                }

                                // assignment_mode = 1, path flow mode
                                {
                                    // assignment_mode = 2, simulated agent flow mode //DTA simulation 

                                    for (int vi = 0; vi < it->second.agent_simu_id_vector.size(); ++vi)
                                    {


                                        int agent_simu_id = it->second.agent_simu_id_vector[vi];
                                        CAgent_Simu* pAgentSimu = g_agent_simu_vector[agent_simu_id];
                                        if (pAgentSimu->agent_id == 81)
                                        {
                                            int idebug = 1;
                                        }
                                        if (assignment.trajectory_diversion_only == 1 && pAgentSimu->diversion_flag == 0)  // diversion flag only, then we skip the non-diversion path
                                            continue;

                                        time_stamp = assignment.g_LoadingStartTimeInMin + pAgentSimu->departure_time_in_min;
                                        for (int nt = 0 + virtual_link_delta; nt < pAgentSimu->path_link_seq_no_vector.size() + 1 - virtual_link_delta; ++nt)
                                        {
                                            float time_in_min = 0;

                                            if (nt < pAgentSimu->path_link_seq_no_vector.size() - virtual_link_delta)
                                                time_in_min = assignment.g_LoadingStartTimeInMin + pAgentSimu->m_Veh_LinkArrivalTime_in_simu_interval[nt] * number_of_seconds_per_interval / 60.0;
                                            else
                                                time_in_min = assignment.g_LoadingStartTimeInMin + pAgentSimu->m_Veh_LinkDepartureTime_in_simu_interval[nt - 1] * number_of_seconds_per_interval / 60.0;  // last link in the path

                                            path_time_vector[nt] = time_in_min;
                                        }

                                        float vehicle_travel_time = pAgentSimu->path_travel_time_in_min;

                                        if (count >= 2000)  // only output up to 2000
                                            break;

                                        // some bugs for output link performances before
                                        fprintf(g_pFilePathMOE, ",%d,%d,%d,%d,%s,%s,%d,%s,1,%.1f,%.4f,%.4f,",
                                            pAgentSimu->agent_id,
                                            g_zone_vector[orig].zone_id,
                                            g_zone_vector[dest].zone_id,
                                            it->second.path_seq_no,
                                            assignment.g_AgentTypeVector[at].display_code.c_str(),
                                            assignment.g_AgentTypeVector[at].agent_type.c_str(),
                                            pAgentSimu->PCE_unit_size,
                                            assignment.g_DemandPeriodVector[tau].demand_period.c_str(),
                                            path_toll,
                                            vehicle_travel_time,
                                            path_distance);

                                        /* Format and print various data */
                                        for (int ni = 0 + virtual_link_delta; ni < pAgentSimu->path_link_seq_no_vector.size(); ++ni)
                                        {
                                            int node_id = g_node_vector[g_link_vector[pAgentSimu->path_link_seq_no_vector[ni]].from_node_seq_no].node_id;
                                            fprintf(g_pFilePathMOE, "%d;", node_id);
                                        }

                                        fprintf(g_pFilePathMOE, ",");

                                        for (int nt = 0 + virtual_link_delta; nt < pAgentSimu->path_link_seq_no_vector.size() + 1 - virtual_link_delta; ++nt)
                                            fprintf(g_pFilePathMOE, "%f;", path_time_vector[nt]);

                                        // time coded
                                        //for (int nt = 0 + virtual_link_delta; nt < pAgentSimu->path_link_seq_no_vector.size() + 1 - virtual_link_delta; ++nt)
                                        //    fprintf(g_pFilePathMOE, "%s;", g_time_coding(path_time_vector[nt]).c_str());

                                        //fprintf(g_pFilePathMOE, "\"LINESTRING (");

                                        //for (int ni = 0 + virtual_link_delta; ni < pAgentSimu->path_link_seq_no_vector.size(); ++ni)
                                        //{

                                        //    int node_no = g_link_vector[pAgentSimu->path_link_seq_no_vector[ni]].from_node_seq_no;
                                        //    fprintf(g_pFilePathMOE, "%f %f", g_node_vector[node_no].x,
                                        //        g_node_vector[node_no].y);

                                        //    if (ni != pAgentSimu->path_link_seq_no_vector.size() - 1)
                                        //        fprintf(g_pFilePathMOE, ", ");
                                        //}


                                        //fprintf(g_pFilePathMOE, ")\"");
                                        fprintf(g_pFilePathMOE, "\n");

                                        count++;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        fclose(g_pFilePathMOE);
    }

    int b_trace_file = false;

    // output trace file
    if (assignment.assignment_mode == 0 || assignment.trajectory_output == 0)  //LUE
    {
        FILE* g_pFilePathMOE = nullptr;
        fopen_ss(&g_pFilePathMOE, "trace.csv", "w");
        fclose(g_pFilePathMOE);
    }
    else if (assignment.assignment_mode >= 1)  //UE mode, or ODME, DTA
    {
        dtalog.output() << "writing trace.csv.." << endl;

        double path_time_vector[MAX_LINK_SIZE_IN_A_PATH];
        FILE* g_pFilePathMOE = nullptr;
        fopen_ss(&g_pFilePathMOE, "trace.csv", "w");

        if (!g_pFilePathMOE)
        {
            dtalog.output() << "File trace.csv cannot be opened." << endl;
            g_program_stop();
        }

        fprintf(g_pFilePathMOE, "agent_id,seq_no,node_id,timestamp,timestamp_in_min,trip_time_in_min,travel_time_in_sec,waiting_time_in_simu_interval,x_coord,y_coord\n");

        int count = 1;

        clock_t start_t, end_t;
        start_t = clock();
        clock_t iteration_t;

        int agent_type_size = assignment.g_AgentTypeVector.size();
        int zone_size = g_zone_vector.size();
        int demand_period_size = assignment.g_DemandPeriodVector.size();

        CColumnVector* p_column_pool;

        float path_toll = 0;
        float path_distance = 0;
        float path_travel_time = 0;
        float time_stamp = 0;

        if (assignment.trajectory_sampling_rate < 0.01)
            assignment.trajectory_sampling_rate = 0.01;
        int sampling_step = 1;

        std::map<int, CColumnPath>::iterator it, it_begin, it_end;

        dtalog.output() << "writing data for " << zone_size << "  zones " << endl;

        for (int orig = 0; orig < zone_size; ++orig)
        {
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
                                if (count % 100000 == 0)
                                {
                                    end_t = clock();
                                    iteration_t = end_t - start_t;
                                    dtalog.output() << "writing " << count / 1000 << "K agents with CPU time " << iteration_t / 1000.0 << " s" << endl;
                                }

                                if (count % sampling_step != 0)
                                    continue;

                                if (count >= 2000)  // only output up to 2000
                                    break;

                                path_toll = 0;
                                path_distance = 0;
                                path_travel_time = 0;
                                path_time_vector[0] = time_stamp;

                                for (int nl = 0; nl < it->second.m_link_size; ++nl)  // arc a
                                {
                                    int link_seq_no = it->second.path_link_vector[nl];
                                    path_toll += g_link_vector[link_seq_no].VDF_period[tau].toll[at];
                                    path_distance += g_link_vector[link_seq_no].link_distance_VDF;
                                    float link_travel_time = g_link_vector[link_seq_no].travel_time_per_period[tau];
                                    path_travel_time += link_travel_time;
                                    time_stamp += link_travel_time;
                                    path_time_vector[nl + 1] = time_stamp;
                                }

                                int virtual_link_delta = 1;

                                // Xuesong: 11/20/2021, need to check again 
                                // fixed routes have physical nodes always, without virtual connectors
                                //if (p_column_pool->bfixed_route)
                                //    virtual_link_delta = 0;

                                // assignment_mode = 1, path flow mode
                                {
                                    // assignment_mode = 2, simulated agent flow mode //DTA simulation 

                                    for (int vi = 0; vi < it->second.agent_simu_id_vector.size(); ++vi)
                                    {

                                        int agent_simu_id = it->second.agent_simu_id_vector[vi];
                                        CAgent_Simu* pAgentSimu = g_agent_simu_vector[agent_simu_id];


                                        if (pAgentSimu->agent_id == 81)
                                        {
                                            int idebug = 1;
                                        }
                                        //if (assignment.trajectory_diversion_only == 1 && pAgentSimu->diversion_flag == 0)  // diversion flag only, then we skip the non-diversion path
                                        //     continue;

                                        for (int nt = 0 + virtual_link_delta; nt < pAgentSimu->path_link_seq_no_vector.size() + 1 - virtual_link_delta; ++nt)
                                        {
                                            float time_in_min = 0;

                                            if (nt < pAgentSimu->path_link_seq_no_vector.size() - virtual_link_delta)
                                                time_in_min = assignment.g_LoadingStartTimeInMin + pAgentSimu->m_Veh_LinkArrivalTime_in_simu_interval[nt] * number_of_seconds_per_interval / 60.0;
                                            else
                                                time_in_min = assignment.g_LoadingStartTimeInMin + pAgentSimu->m_Veh_LinkDepartureTime_in_simu_interval[nt - 1] * number_of_seconds_per_interval / 60.0;  // last link in the path

                                            path_time_vector[nt] = time_in_min;
                                        }


                                        for (int nt = 0 + virtual_link_delta; nt < pAgentSimu->path_link_seq_no_vector.size() - virtual_link_delta; ++nt)
                                        {
                                            float trip_time_in_min = path_time_vector[nt] - path_time_vector[0 + virtual_link_delta];

                                            int node_id = g_node_vector[g_link_vector[pAgentSimu->path_link_seq_no_vector[nt]].from_node_seq_no].node_id;

                                            if (nt >= pAgentSimu->path_link_seq_no_vector.size() - virtual_link_delta)
                                            {
                                                node_id = g_node_vector[g_link_vector[pAgentSimu->path_link_seq_no_vector[nt]].to_node_seq_no].node_id;

                                            }
                                            int link_seq_no = it->second.path_link_vector[nt];
                                            float travel_time_in_sec = (path_time_vector[nt + 1] - path_time_vector[nt]) * 60;
                                            int waiting_time_in_simu_interaval = (path_time_vector[nt + 1] - path_time_vector[nt] - g_link_vector[link_seq_no].free_flow_travel_time_in_min) * number_of_simu_intervals_in_min;

                                            int node_no = g_link_vector[pAgentSimu->path_link_seq_no_vector[nt]].from_node_seq_no;

                                            fprintf(g_pFilePathMOE, "%d,%d,%d,T%s,%.5f,%f,%.4f,%d,%f,%f\n",
                                                pAgentSimu->agent_id, nt, node_id,
                                                g_time_coding(path_time_vector[nt]).c_str(),
                                                path_time_vector[nt],
                                                trip_time_in_min,
                                                travel_time_in_sec,
                                                waiting_time_in_simu_interaval,
                                                g_node_vector[node_no].x, g_node_vector[node_no].y);

                                        }
                                        count++;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        fclose(g_pFilePathMOE);
    }


}

void g_output_trajectory_bin(Assignment& assignment)
{

    if (assignment.assignment_mode == 0 || assignment.trajectory_output == 0)  //LUE
    {
        FILE* g_pFilePathMOE = nullptr;
        fopen_ss(&g_pFilePathMOE, "trajectory.bin", "wb");
        fclose(g_pFilePathMOE);
    }
    else if (assignment.assignment_mode >= 1)  //UE mode, or ODME, DTA
    {
        dtalog.output() << "writing trajectory.bin.." << endl;

        int path_node_vector[MAX_LINK_SIZE_IN_A_PATH];
        double path_time_vector[MAX_LINK_SIZE_IN_A_PATH];
        FILE* g_pFilePathMOE = nullptr;
        fopen_ss(&g_pFilePathMOE, "trajectory.bin", "wb");

        if (!g_pFilePathMOE)
        {
            dtalog.output() << "File trajectory.bin cannot be opened." << endl;
            g_program_stop();
        }


        struct STrajectoryHeader
        {
            int agent_id, o_zone_id, d_zone_id, path_id, display_code, agent_type, PCE_unit, demand_period;
            int node_size, link_size;
            float volume, toll, travel_time, distance_km;
        };

        STrajectoryHeader header;

        int count = 1;

        clock_t start_t, end_t;
        start_t = clock();
        clock_t iteration_t;

        int agent_type_size = assignment.g_AgentTypeVector.size();
        int zone_size = g_zone_vector.size();
        int demand_period_size = assignment.g_DemandPeriodVector.size();

        CColumnVector* p_column_pool;

        float path_toll = 0;
        float path_distance = 0;
        float path_travel_time = 0;
        float time_stamp = 0;

        if (assignment.trajectory_sampling_rate < 0.01)
            assignment.trajectory_sampling_rate = 0.01;
        int sampling_step = 100 / int(100 * assignment.trajectory_sampling_rate + 0.5);

        std::map<int, CColumnPath>::iterator it, it_begin, it_end;

        dtalog.output() << "writing data for " << zone_size << "  zones " << endl;

        for (int orig = 0; orig < zone_size; ++orig)
        {
            if (g_zone_vector[orig].zone_id % 100 == 0)
                dtalog.output() << "o zone id =  " << g_zone_vector[orig].zone_id << endl;

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
                                if (count % 100000 == 0)
                                {
                                    end_t = clock();
                                    iteration_t = end_t - start_t;
                                    dtalog.output() << "writing " << count / 1000 << "K binary agents with CPU time " << iteration_t / 1000.0 << " s" << endl;
                                }

                                if (count % sampling_step != 0)
                                    continue;


                                path_toll = 0;
                                path_distance = 0;
                                path_travel_time = 0;
                                path_time_vector[0] = time_stamp;

                                for (int nl = 0; nl < it->second.m_link_size; ++nl)  // arc a
                                {
                                    int link_seq_no = it->second.path_link_vector[nl];
                                    path_toll += g_link_vector[link_seq_no].VDF_period[tau].toll[at];
                                    path_distance += g_link_vector[link_seq_no].link_distance_VDF;
                                    float link_travel_time = g_link_vector[link_seq_no].travel_time_per_period[tau];

                                    time_stamp += link_travel_time;
                                    path_time_vector[nl + 1] = time_stamp;
                                }

                                int virtual_link_delta = 1;
                                int virtual_begin_link_delta = 1;
                                int virtual_end_link_delta = 1;
                                // fixed routes have physical nodes always, without virtual connectors
                                if (p_column_pool->bfixed_route)
                                {
                                    virtual_begin_link_delta = 0;
                                    virtual_end_link_delta = 1;

                                }

                                // assignment_mode = 1, path flow mode
                                {
                                    // assignment_mode = 2, simulated agent flow mode //DTA simulation 

                                    for (int vi = 0; vi < it->second.agent_simu_id_vector.size(); ++vi)
                                    {


                                        int agent_simu_id = it->second.agent_simu_id_vector[vi];
                                        CAgent_Simu* pAgentSimu = g_agent_simu_vector[agent_simu_id];
                                        if (pAgentSimu->agent_id == 81)
                                        {
                                            int idebug = 1;
                                        }
                                        if (assignment.trajectory_diversion_only == 1 && pAgentSimu->diversion_flag == 0)  // diversion flag only, then we skip the non-diversion path
                                            continue;

                                        time_stamp = assignment.g_LoadingStartTimeInMin + pAgentSimu->departure_time_in_min;
                                        for (int nt = 0 + virtual_link_delta; nt < pAgentSimu->path_link_seq_no_vector.size() + 1 - virtual_link_delta; ++nt)
                                        {
                                            double time_in_min = 0;

                                            if (nt < pAgentSimu->path_link_seq_no_vector.size() - virtual_link_delta)
                                                time_in_min = assignment.g_LoadingStartTimeInMin + pAgentSimu->m_Veh_LinkArrivalTime_in_simu_interval[nt] * number_of_seconds_per_interval / 60.0;
                                            else
                                                time_in_min = assignment.g_LoadingStartTimeInMin + pAgentSimu->m_Veh_LinkDepartureTime_in_simu_interval[nt - 1] * number_of_seconds_per_interval / 60.0;  // last link in the path

                                            path_time_vector[nt- virtual_link_delta] = time_in_min;
                                        }

                                        float vehicle_travel_time = pAgentSimu->path_travel_time_in_min;


                                        header.agent_id = pAgentSimu->agent_id;
                                        header.o_zone_id = g_zone_vector[orig].zone_id;
                                        header.d_zone_id = g_zone_vector[dest].zone_id;
                                        header.path_id = it->second.path_seq_no;
                                        header.display_code = 0; // assignment.g_AgentTypeVector[at].display_code_no;
                                        header.agent_type = 0; // assignment.g_AgentTypeVector[at].agent_type_no;
                                        header.PCE_unit = 0; // pAgentSimu->PCE_unit_size;
                                        header.demand_period = tau;
                                        header.toll = path_toll;
                                        header.travel_time = vehicle_travel_time;
                                        header.distance_km = path_distance;
                                        header.node_size = pAgentSimu->path_link_seq_no_vector.size() - virtual_link_delta;
                                        header.link_size = pAgentSimu->path_link_seq_no_vector.size() + 1 - virtual_link_delta;

                                        fwrite(&header, sizeof(struct STrajectoryHeader), 1, g_pFilePathMOE);

                    /* Format and print various data */
                              
                                            for (int ni = 0 + virtual_link_delta; ni < pAgentSimu->path_link_seq_no_vector.size(); ++ni)
                                            {
                                                path_node_vector[ni- virtual_link_delta] = g_node_vector[g_link_vector[pAgentSimu->path_link_seq_no_vector[ni]].from_node_seq_no].node_id;
                                            }

                                            fwrite(&path_node_vector, sizeof(int), header.node_size, g_pFilePathMOE);
                                            fwrite(&path_time_vector, sizeof(double), header.link_size, g_pFilePathMOE);


                                        count++;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

        end_t = clock();
        iteration_t = end_t - start_t;
        dtalog.output() << "Comlete writing " << count / 1000 << "K binary agents with CPU time " << iteration_t / 1000.0 << " s." << endl;
        fclose(g_pFilePathMOE);
    }

}

void g_output_simulation_result(Assignment& assignment)
{
    g_output_dynamic_link_performance(assignment, 1);
    g_output_dynamic_link_performance(assignment, 2);


    if (assignment.assignment_mode == 2)  //DTA mode
    {
        g_output_dynamic_link_state(assignment, 1);
    }

    g_output_trajectory_bin(assignment);

    //g_output_simulation_agents(assignment);
    g_output_trajectory_csv(assignment);
}

void g_OutputModelFiles(int mode)
{
    if (mode == 1)
    {
        FILE* g_pFileModelNode = fopen("model_node.csv", "w");

        if (g_pFileModelNode != NULL)
        {
            fprintf(g_pFileModelNode, "node_id,node_no,node_type,#_of_outgoing_nodes,activity_node_flag,agent_type,zone_id,cell_id,cell_code,info_zone_flag,x_coord,y_coord\n");
            for (int i = 0; i < g_node_vector.size(); i++)
            {


                if (g_node_vector[i].node_id >= 0)  //for all physical links
                {

                    fprintf(g_pFileModelNode, "%d,%d,%s,%d,%d,%s,%d,%ld,%s,%d,%f,%f\n",
                        g_node_vector[i].node_id,
                        g_node_vector[i].node_seq_no,
                        g_node_vector[i].node_type.c_str(),
                        g_node_vector[i].m_outgoing_link_seq_no_vector.size(),
                        g_node_vector[i].is_activity_node,
                        g_node_vector[i].agent_type_str.c_str(),
                        g_node_vector[i].zone_org_id,
                        g_node_vector[i].cell_id,
                        g_node_vector[i].cell_str.c_str(),
                        g_node_vector[i].is_information_zone,
                        g_node_vector[i].x,
                        g_node_vector[i].y
                    );

                }

            }

            fclose(g_pFileModelNode);
        }
        else
        {
            dtalog.output() << "Error: File model_node.csv cannot be opened.\n It might be currently used and locked by EXCEL." << endl;
            g_program_stop();


        }

    }

    if (mode == 2)
    {
        FILE* g_pFileModelLink = fopen("model_link.csv", "w");

        if (g_pFileModelLink != NULL)
        {
            fprintf(g_pFileModelLink, "link_id,link_no,from_node_id,to_node_id,link_type,link_type_name,lanes,link_distance_VDF,free_speed,fftt,capacity,allow_uses,geometry\n");

            //VDF_fftt1,VDF_cap1,VDF_alpha1,VDF_beta1
            for (int i = 0; i < g_link_vector.size(); i++)
            {

                fprintf(g_pFileModelLink, "%s,%d,%d,%d,%d,%s,%d,%f,%f,%f,%f,%s,",
                    g_link_vector[i].link_id.c_str(),
                    g_link_vector[i].link_seq_no,
                    g_node_vector[g_link_vector[i].from_node_seq_no].node_id,
                    g_node_vector[g_link_vector[i].to_node_seq_no].node_id,
                    g_link_vector[i].link_type,
                    g_link_vector[i].link_type_name.c_str(),
                    g_link_vector[i].number_of_lanes,
                    g_link_vector[i].link_distance_VDF,
                    g_link_vector[i].free_speed,
                    g_link_vector[i].free_flow_travel_time_in_min,
                    g_link_vector[i].lane_capacity,
                    g_link_vector[i].VDF_period[0].allowed_uses.c_str()
                    //g_link_vector[i].VDF_period[0].FFTT,
                   //g_link_vector[i].VDF_period[0].period_capacity,
                   //g_link_vector[i].VDF_period[0].alpha,
                   //g_link_vector[i].VDF_period[0].beta,
                );

                if (g_link_vector[i].geometry.size() > 0)
                {
                    fprintf(g_pFileModelLink, "\"%s\",", g_link_vector[i].geometry.c_str());
                }
                else
                {
                    fprintf(g_pFileModelLink, "\"LINESTRING (");

                    fprintf(g_pFileModelLink, "%f %f,", g_node_vector[g_link_vector[i].from_node_seq_no].x, g_node_vector[g_link_vector[i].from_node_seq_no].y);
                    fprintf(g_pFileModelLink, "%f %f", g_node_vector[g_link_vector[i].to_node_seq_no].x, g_node_vector[g_link_vector[i].to_node_seq_no].y);

                    fprintf(g_pFileModelLink, ")\"");
                }

                fprintf(g_pFileModelLink, "\n");

            }


            fclose(g_pFileModelLink);
        }
        else
        {
            dtalog.output() << "Error: File model_link.csv cannot be opened.\n It might be currently used and locked by EXCEL." << endl;
            g_program_stop();

        }

    }
    if (mode == 3)
    {
        int connector_count = 0;
        for (int i = 0; i < g_link_vector.size(); i++)
        {
            if (g_link_vector[i].link_type == 1000)  // connector
            {
                connector_count += 1;
            }
        }

        if(connector_count >=1)
        {
        FILE* g_pFileModelLink = fopen("access_link.csv", "w");

        if (g_pFileModelLink != NULL)
        {
            fprintf(g_pFileModelLink, "link_id,link_no,from_node_id,to_node_id,link_type,link_type_name,lanes,link_distance_VDF,free_speed,fftt,capacity,allow_uses,geometry\n");

            //VDF_fftt1,VDF_cap1,VDF_alpha1,VDF_beta1
            for (int i = 0; i < g_link_vector.size(); i++)
            {
                if (g_link_vector[i].link_type == 1000 )  // connector
                {

                    fprintf(g_pFileModelLink, "%s,%d,%d,%d,%d,%s,%d,%f,%f,%f,%f,%s,",
                        g_link_vector[i].link_id.c_str(),
                        g_link_vector[i].link_seq_no,
                        g_node_vector[g_link_vector[i].from_node_seq_no].node_id,
                        g_node_vector[g_link_vector[i].to_node_seq_no].node_id,
                        g_link_vector[i].link_type,
                        g_link_vector[i].link_type_name.c_str(),
                        g_link_vector[i].number_of_lanes,
                        g_link_vector[i].link_distance_VDF,
                        g_link_vector[i].free_speed,
                        g_link_vector[i].free_flow_travel_time_in_min,
                        g_link_vector[i].lane_capacity,
                        g_link_vector[i].VDF_period[0].allowed_uses.c_str()
                         //g_link_vector[i].VDF_period[0].FFTT,
                        //g_link_vector[i].VDF_period[0].period_capacity,
                        //g_link_vector[i].VDF_period[0].alpha,
                        //g_link_vector[i].VDF_period[0].beta,
                    );

                    if (g_link_vector[i].geometry.size() > 0)
                    {
                        fprintf(g_pFileModelLink, "\"%s\",", g_link_vector[i].geometry.c_str());
                    }else
                    {
                    fprintf(g_pFileModelLink, "\"LINESTRING (");

                    fprintf(g_pFileModelLink, "%f %f,", g_node_vector[g_link_vector[i].from_node_seq_no].x, g_node_vector[g_link_vector[i].from_node_seq_no].y);
                    fprintf(g_pFileModelLink, "%f %f", g_node_vector[g_link_vector[i].to_node_seq_no].x, g_node_vector[g_link_vector[i].to_node_seq_no].y);

                    fprintf(g_pFileModelLink, ")\"");
                    }

                    fprintf(g_pFileModelLink, "\n");

                }


            }

            fclose(g_pFileModelLink);
        }
        else
        {
            dtalog.output() << "Error: File access_link.csv cannot be opened.\n It might be currently used and locked by EXCEL." << endl;
            g_program_stop();

        }

        }

    }

    if (mode == 3)  // cell
    {
        //FILE* g_pFileZone = nullptr;
        //g_pFileZone = fopen("model_cell.csv", "w");

        //if (g_pFileZone == NULL)
        //{
        //    cout << "File model_cell.csv cannot be opened." << endl;
        //    g_program_stop();
        //}
        //else
        //{


        //    fprintf(g_pFileZone, "cell_code,geometry\n");

        //    std::map<string, CInfoCell>::iterator it;

        //    for (it = g_info_cell_map.begin(); it != g_info_cell_map.end(); ++it)
        //    {

        //        fprintf(g_pFileZone, "%s,", it->first.c_str());
        //        fprintf(g_pFileZone, "\"LINESTRING (");

        //        for (int s = 0; s < it->second.m_ShapePoints.size(); s++)
        //        {
        //            fprintf(g_pFileZone, "%f %f,", it->second.m_ShapePoints[s].x, it->second.m_ShapePoints[s].y);
        //        }

        //        fprintf(g_pFileZone, ")\"");
        //        fprintf(g_pFileZone, "\n");
        //    }
        //    fclose(g_pFileZone);
        //}
    }

    if (mode == 10)
    {
        FILE* g_pFileModel_LC = fopen("model_shortest_path_tree.csv", "w");

        if (g_pFileModel_LC != NULL)
        {
            fprintf(g_pFileModel_LC, "iteration,agent_type,zone_id,node_id,pred,label_cost\n");
            for (int i = 0; i < g_node_vector.size(); i++)
            {


                if (g_node_vector[i].node_id >= 0)  //for all physical links
                {
                    std::map<string, float> ::iterator it;

                    for (it = g_node_vector[i].label_cost_per_iteration_map.begin(); it != g_node_vector[i].label_cost_per_iteration_map.end(); ++it)
                    {
                        int node_pred_id = -1;
                        int pred_no = g_node_vector[i].pred_per_iteration_map[it->first];
                        if (pred_no >= 0)
                            node_pred_id = g_node_vector[pred_no].node_id;

                        //if(it->second < 100000)
                        {
                        fprintf(g_pFileModel_LC, "%s,%d,%f,\n", it->first.c_str(), node_pred_id, it->second);
                        }
                    }

                }

            }

            fclose(g_pFileModel_LC);
        }
        else
        {
            dtalog.output() << "Error: File model_label_cost_tree.csv cannot be opened.\n It might be currently used and locked by EXCEL." << endl;
            g_program_stop();


        }

    }

}



// FILE* g_pFileOutputLog = nullptr;

