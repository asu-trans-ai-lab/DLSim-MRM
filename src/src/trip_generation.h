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

double g_random_generate_activity_nodes(Assignment& assignment)
{

    int activity_node_count = 0;                    // check how many actiity nodes
    for (int i = 0; i < g_node_vector.size(); i++)  
    {

        if (g_node_vector[i].is_activity_node >= 1)
        {
            activity_node_count++;
        }
    }


    if (activity_node_count <= 1    )       // random generation of activity locations 
    {
        activity_node_count = 0;
        int sampling_rate = 10;

        for (int i = 0; i < g_node_vector.size(); i++)
        {

            if (i % sampling_rate == 0)
            {
                g_node_vector[i].is_activity_node = 10;//random generation
                activity_node_count++;
            }
        }

        if (activity_node_count <= 1)
        {
            activity_node_count = 0;
            sampling_rate = 2;

            for (int i = 0; i < g_node_vector.size(); i++)
            {

                if (i % sampling_rate == 0)
                {
                    g_node_vector[i].is_activity_node = 10;//random generation
                    activity_node_count++;
                }
            }
            // still no activity nodes, define all nodes as activity nodes
            if (activity_node_count <= 1)
            {
                activity_node_count = 0;

                for (int i = 0; i < g_node_vector.size(); i++)
                {

                    g_node_vector[i].is_activity_node = 10;//random generation
                    activity_node_count++;
                }
            }
        }

    }


    // calculate avg near by distance; 
    double total_near_by_distance = 0;
    activity_node_count = 0;
    for (int i = 0; i < g_node_vector.size(); i++)
    {
        double min_near_by_distance = 100;
        if (g_node_vector[i].is_activity_node)
        {
            activity_node_count++;
            for (int j = 0; j < g_node_vector.size(); j++)
            {
                if (i != j && g_node_vector[j].is_activity_node)
                {
                    double near_by_distance = g_calculate_p2p_distance_in_mile_from_latitude_longitude(g_node_vector[i].x, g_node_vector[i].y, g_node_vector[j].x, g_node_vector[j].y);

                    if (near_by_distance < min_near_by_distance)
                        min_near_by_distance = near_by_distance;

                }

            }

            total_near_by_distance += min_near_by_distance;
            activity_node_count++;
        }
    }

    double nearby_distance = total_near_by_distance / max(1, activity_node_count);
    return nearby_distance;

}

void g_grid_zone_generation(Assignment& assignment)
{
    dtalog.output() << "Step 1.4.1: QEM mode for creating node 2 zone mapping" << endl;

    double activity_nearbydistance = g_random_generate_activity_nodes(assignment);
    // initialization of grid rectangle boundary
    double left = 100000000;
    double right = -100000000;
    double top = -1000000000;
    double  bottom = 1000000000;

    for (int i = 0; i < g_node_vector.size(); i++)
    {
        // exapnd the grid boundary according to the nodes
        left = min(left, g_node_vector[i].x);
        right = max(right, g_node_vector[i].x);
        top = max(top, g_node_vector[i].y);
        bottom = min(bottom, g_node_vector[i].y);

    }

    int grid_size = 8;

    if (g_node_vector.size() > 3000)
        grid_size = 10;
    if (g_node_vector.size() > 10000)
        grid_size = 20;
    if (g_node_vector.size() > 40000)
        grid_size = 30;

    double temp_resolution = (((right - left) / grid_size + (top - bottom) / grid_size)) / 2.0;

    if (activity_nearbydistance * 4 < temp_resolution)
    {
        temp_resolution = activity_nearbydistance * 4;

    }


    vector<double> ResolutionVector;

    ResolutionVector.push_back(0.00005);
    ResolutionVector.push_back(0.0001);
    ResolutionVector.push_back(0.0002);
    ResolutionVector.push_back(0.00025);
    ResolutionVector.push_back(0.0005);
    ResolutionVector.push_back(0.00075);
    ResolutionVector.push_back(0.001);
    ResolutionVector.push_back(0.002);
    ResolutionVector.push_back(0.0025);
    ResolutionVector.push_back(0.005);
    ResolutionVector.push_back(0.0075);
    ResolutionVector.push_back(0.01);
    ResolutionVector.push_back(0.02);
    ResolutionVector.push_back(0.025);
    ResolutionVector.push_back(0.05);
    ResolutionVector.push_back(0.075);
    ResolutionVector.push_back(0.1);
    ResolutionVector.push_back(0.2);
    ResolutionVector.push_back(0.25);
    ResolutionVector.push_back(0.5);
    ResolutionVector.push_back(0.75);
    ResolutionVector.push_back(1);
    ResolutionVector.push_back(2);
    ResolutionVector.push_back(2.5);
    ResolutionVector.push_back(5);
    ResolutionVector.push_back(7.5);
    ResolutionVector.push_back(10);
    ResolutionVector.push_back(20);
    ResolutionVector.push_back(25);
    ResolutionVector.push_back(50);
    ResolutionVector.push_back(75);

    double ClosestResolution = 1;

    if (temp_resolution < ResolutionVector[0])
        temp_resolution = ResolutionVector[0];

    for (unsigned int i = 0; i < ResolutionVector.size() - 1; i++)
    {
        if ((temp_resolution > ResolutionVector[i] + 0.000001) && temp_resolution < ResolutionVector[i + 1])
        {
            temp_resolution = ResolutionVector[i + 1]; // round up
            break;

        }
    }

    assignment.m_GridResolution = temp_resolution;

    assignment.zone_id_2_node_no_mapping.clear();
    dtalog.output() << "Step 1.4.2: Grid Resolution " << assignment.m_GridResolution << endl;

    int activity_node_count = 0;
    for (int i = 0; i < g_node_vector.size(); i++)
    {

        if (g_node_vector[i].is_activity_node >= 1)
        {

            if (g_node_vector[i].node_id == 966)
            {
                int itest = 1;
            }
            __int64 cell_id = g_GetCellID(g_node_vector[i].x, g_node_vector[i].y, assignment.m_GridResolution);
            int zone_id;

            if (assignment.cell_id_mapping.find(cell_id) == assignment.cell_id_mapping.end())  // create a cell
            {
                //create zone
                assignment.cell_id_mapping[cell_id] = g_node_vector[i].node_id;


                dtalog.output() << "Step 1.2: creating cell " << cell_id << " using node id " << g_node_vector[i].node_id << endl;

                zone_id = assignment.cell_id_mapping[cell_id]; // which is the node id when a cell is created. 
                if (assignment.zone_id_2_node_no_mapping.find(zone_id) == assignment.zone_id_2_node_no_mapping.end()) // create a zone 
                {
                    dtalog.output() << "Step 1.2: creating zone " << zone_id << " using node id " << g_node_vector[i].node_id << endl;
                    //create zone
                    assignment.zone_id_2_node_no_mapping[zone_id] = i;
                    assignment.zone_id_2_cell_id_mapping[zone_id] = cell_id;
                    g_node_vector[i].zone_org_id = zone_id;
                }
            }
            else
            {
                zone_id = assignment.cell_id_mapping[cell_id]; // which is the node id when a cell is created. 
                // for physcial nodes because only centriod can have valid zone_id.
                g_node_vector[i].zone_org_id = zone_id;

            }

            activity_node_count++;


        }
    }

    dtalog.output() << "Step 1.4.3: creating " << assignment.zone_id_2_node_no_mapping.size() << " zones." << " # of activity nodes =" << activity_node_count << endl;

}

void g_create_zone_vector(Assignment& assignment)
{
    std::map<int, int> waring_message_link_type_map;
    // initialize zone vector
    dtalog.output() << "Step 1.5: Initializing O-D zone vector..." << endl;

    std::map<int, int>::iterator it;

    for (it = assignment.zone_id_2_node_no_mapping.begin(); it != assignment.zone_id_2_node_no_mapping.end(); ++it)
    {
        COZone ozone;

        if (it->first == 966)
        {
            int itest = 1;
        }
        // for each zone, we have to also create centriod
        ozone.zone_id = it->first;  // zone_id
        ozone.cell_id = assignment.zone_id_2_cell_id_mapping[it->first];
        ozone.zone_seq_no = g_zone_vector.size();
        ozone.cell_x = g_node_vector[it->second].x;
        ozone.cell_y = g_node_vector[it->second].y;

        dtalog.output() << "create zone id = " << ozone.zone_id << " with representive node id " << it->second << ",x = " << g_node_vector[it->second].x << ",y=" <<
            ozone.cell_y << endl;


        assignment.g_zoneid_to_zone_seq_no_mapping[ozone.zone_id] = ozone.zone_seq_no;  // create the zone id to zone seq no mapping

        // create a centriod
        CNode node;
        // very large number as a special id
        node.node_id = -1 * (ozone.zone_id) - 1000000;
        node.node_seq_no = g_node_vector.size();
        assignment.g_node_id_to_seq_no_map[node.node_id] = node.node_seq_no;
        node.zone_id = ozone.zone_id;
        // push it to the global node vector
        g_node_vector.push_back(node);
        assignment.g_number_of_nodes++;

        ozone.node_seq_no = node.node_seq_no;
        // this should be the only one place that defines this mapping
        assignment.zone_id_to_centriod_node_no_mapping[ozone.zone_id] = node.node_seq_no;
        // add element into vector
        g_zone_vector.push_back(ozone);
    }

}
void g_trip_generation(Assignment& assignment)
{
    // accessibility 
    for (int tau = 0; tau < assignment.g_DemandPeriodVector.size(); ++tau)
    {
        for (int at = 0; at < assignment.g_AgentTypeVector.size(); ++at)
        {

            for (int orig = 0; orig < g_zone_vector.size(); ++orig)  // o
            {
                for (int dest = 0; dest < g_zone_vector.size(); ++dest)  // d
                {
                    if (orig != dest)
                    {
                        float distance_in_mile = g_calculate_p2p_distance_in_mile_from_latitude_longitude(g_zone_vector[orig].cell_x, g_zone_vector[orig].cell_y, g_zone_vector[dest].cell_x, g_zone_vector[dest].cell_y);
                        g_zone_vector[orig].m_ODAccessibilityMatrix[at][tau].distance_map[dest] = distance_in_mile;
                        float travel_time_in_min = distance_in_mile/ 30 * 60;  // default speed as 30 miles per hour

                        g_zone_vector[orig].m_ODAccessibilityMatrix[at][tau].value_map[dest] = travel_time_in_min;
                    }

                }

            }
        }
    }

    int out_of_bound_log_count = 0;
    int trip_accessibility_log_count = 0;
    int trip_distribution_log_count = 0;

    for (int tau = 0; tau < assignment.g_DemandPeriodVector.size(); ++tau)
    {
        for (int at = 0; at < assignment.g_AgentTypeVector.size(); ++at)
        {
            for (int orig = 0; orig < g_zone_vector.size(); ++orig)  // o
                                                                     // gravity model;
            {
                if (g_zone_vector[orig].gravity_production[at] > 0.00001)
                {

                    float total_attraction_utility = 0;
                    int count = 0;

                    for (int d = 0; d < g_zone_vector.size(); ++d)
                    {
                        if (orig != d)
                        {

                            g_zone_vector[d].gravity_est_attraction[at] = 0;

//                            float cut_off = assignment.g_AgentTypeVector[at].trip_time_budget_in_min;
                            float cut_off = 40;
                            if (g_zone_vector[d].gravity_attraction[at] > 0)
                            {
                                //double disutility = g_zone_vector[orig].m_ODAccessibilityMatrix[at][tau].value_map[d] * beta;
                                double exp_disutility = g_zone_vector[d].gravity_attraction[at];
                                g_zone_vector[orig].m_ODAccessibilityMatrix[at][tau].disutility_map[d] = exp_disutility;

                                if (g_zone_vector[orig].m_ODAccessibilityMatrix[at][tau].value_map[d] < cut_off)
                                {
                                    if (trip_accessibility_log_count <= 100)
                                    {
                                        dtalog.output() << "agent type: " << assignment.g_AgentTypeVector[at].agent_type.c_str() << ", o: " << orig << ",d:" << d <<
                                            ", gc distance = " << g_zone_vector[orig].m_ODAccessibilityMatrix[at][tau].distance_map[d] <<
                                            ", travel time =" << g_zone_vector[orig].m_ODAccessibilityMatrix[at][tau].value_map[d] <<
                                            ",value = " << exp_disutility << endl;
                                    }
                                    total_attraction_utility += exp_disutility;
                                    trip_accessibility_log_count++;
                                    count++;
                                }
                                else
                                {
                                    if (out_of_bound_log_count < 10)
                                    {
                                        dtalog.output() << "out of bound: agent type: " << assignment.g_AgentTypeVector[at].agent_type.c_str() << ",o: " << orig << ",d:" << d <<
                                            ", gc distance = " << g_zone_vector[orig].m_ODAccessibilityMatrix[at][tau].distance_map[d] <<
                                            ", travel time =" << g_zone_vector[orig].m_ODAccessibilityMatrix[at][tau].value_map[d] <<
                                            ",value = " << exp_disutility << endl;
                                    }

                                    out_of_bound_log_count++;

                                }
                            }
                        }
                    }

                    dtalog.output() << "o: " << orig << ", total_attraction_utility =" << total_attraction_utility << endl;

                    if (count > 0)
                    {
                        for (int dest = 0; dest < g_zone_vector.size(); ++dest)  // d
                        {
                            if (orig != dest)
                            {
                                if (g_zone_vector[dest].gravity_attraction[at] > 0)
                                {
//                                    float cut_off = assignment.g_AgentTypeVector[at].trip_time_budget_in_min;
                                    float cut_off = 40;
                                    if (g_zone_vector[orig].m_ODAccessibilityMatrix[at][tau].value_map[dest] < cut_off)
                                    {

                                        double ratio = g_zone_vector[orig].m_ODAccessibilityMatrix[at][tau].disutility_map[dest] / total_attraction_utility;


                                        g_zone_vector[orig].m_ODMatrix[at][tau].value_map[dest] = g_zone_vector[orig].gravity_production[at] * ratio;

                                        if (trip_distribution_log_count < 100)
                                        {
                                            dtalog.output() << "agent type: " << assignment.g_AgentTypeVector[at].agent_type.c_str() << ", o: " << orig << ",d:" << dest << ", ratio =" << ratio <<
                                                ",trip = " << g_zone_vector[orig].m_ODMatrix[at][tau].value_map[dest] << endl;
                                        }
                                        trip_distribution_log_count++;

                                        g_zone_vector[dest].gravity_est_attraction[at] += g_zone_vector[orig].m_ODMatrix[at][tau].value_map[dest];
                                    }
                                }

                            }
                        }
                    }
                }

            }
        }
    }
}
void g_writing_demand_files(Assignment& assignment)
{
dtalog.output() << "writing demand_geo.csv.." << endl;

FILE* g_pFileODMatrix = nullptr;
fopen_ss(&g_pFileODMatrix, "demand_geo.csv", "w");

if (!g_pFileODMatrix)
{
    dtalog.output() << "File demand_geo.csv cannot be opened." << endl;
    g_program_stop();
}
else
{

    fprintf(g_pFileODMatrix, "demand_period,time_period,agent_type,o_zone_id,d_zone_id,volume,geometry\n");
    int demand_writing_log_count = 0;
    // reset the estimated production and attraction
    for (int orig = 0; orig < g_zone_vector.size(); ++orig)  // o
    {
        for (int dest = 0; dest < g_zone_vector.size(); ++dest)  // d
        {
            for (int tau = 0; tau < assignment.g_DemandPeriodVector.size(); ++tau)
            {
                for (int at = 0; at < assignment.g_AgentTypeVector.size(); ++at)
                {
                    if (g_zone_vector[orig].gravity_production[at] >= 0)
                    {
                        if (g_zone_vector[dest].gravity_attraction[at] > 0)
                        {
                            float value = 0;
                            if (g_zone_vector[orig].m_ODMatrix[at][tau].value_map.find(dest) != g_zone_vector[orig].m_ODMatrix[at][tau].value_map.end())
                            {
                                value = g_zone_vector[orig].m_ODMatrix[at][tau].value_map[dest];
                            }

                            if (value > 0.000001)
                            {
                                fprintf(g_pFileODMatrix, "%s,%s,%s,%d,%d,%.4f,", assignment.g_DemandPeriodVector[tau].demand_period.c_str(), assignment.g_DemandPeriodVector[tau].time_period.c_str(),

                                    assignment.g_AgentTypeVector[at].agent_type.c_str(), g_zone_vector[orig].zone_id, g_zone_vector[dest].zone_id, value);

                                if (demand_writing_log_count < 100)
                                {
                                    dtalog.output() << "orig= " << g_zone_vector[orig].zone_id << " dest= " << g_zone_vector[dest].zone_id << ":" << value << endl;
                                }
                                demand_writing_log_count++;

                                fprintf(g_pFileODMatrix, "\"LINESTRING (");

                                fprintf(g_pFileODMatrix, "%f %f,", g_zone_vector[orig].cell_x, g_zone_vector[orig].cell_y);
                                fprintf(g_pFileODMatrix, "%f %f,", g_zone_vector[dest].cell_x, g_zone_vector[dest].cell_y);
                                fprintf(g_pFileODMatrix, ")\"");
                                fprintf(g_pFileODMatrix, "\n");
                            }

                        }
                    }

                }
            }
        }
    }

    fclose(g_pFileODMatrix);
}


////////////////////////////
dtalog.output() << "writing input_matrix.csv.." << endl;


fopen_ss(&g_pFileODMatrix, "input_matrix.csv", "w");

if (!g_pFileODMatrix)
{
    dtalog.output() << "File input_matrix.csv cannot be opened." << endl;
    g_program_stop();
}
else
{
    for (int tau = 0; tau < assignment.g_DemandPeriodVector.size(); ++tau)
    {
        for (int at = 0; at < assignment.g_AgentTypeVector.size(); ++at)
        {

            fprintf(g_pFileODMatrix, "od,");  // first line
            for (int d = 0; d < g_zone_vector.size(); ++d)
            {
                fprintf(g_pFileODMatrix, "%d,", g_zone_vector[d].zone_id);
            }

            fprintf(g_pFileODMatrix, "\n");  // return key at the first line

//            fprintf(g_pFileODMatrix, "subtotal_est,subtotal_target,ratio_diff\n");

            // reset the estimated production and attraction
            for (int orig = 0; orig < g_zone_vector.size(); ++orig)  // o
            {
////                fprintf(g_pFileODMatrix, "%s,%s,", assignment.g_DemandPeriodVector[tau].demand_period.c_str(),
//                    assignment.g_AgentTypeVector[at].agent_type.c_str());
                float total_production = 0;

                fprintf(g_pFileODMatrix, "%d,", g_zone_vector[orig].zone_id);  // origin zone id 

                for (int dest = 0; dest < g_zone_vector.size(); ++dest)  // d
                {
                    float value = 0;
                    if (g_zone_vector[orig].gravity_production[at] >= 0 && g_zone_vector[dest].gravity_attraction[at] > 0)
                    {
                        value = g_zone_vector[orig].m_ODMatrix[at][tau].value_map[dest];
                    }
                    total_production += value;
                    fprintf(g_pFileODMatrix, "%f,", value);

                }
                fprintf(g_pFileODMatrix, "\n");  // return key at each line 
                //                float percentage_difference = (total_production - g_zone_vector[orig].gravity_production[at]) / max(0.001, g_zone_vector[orig].gravity_production[at]);
//                fprintf(g_pFileODMatrix, "%f,%f,%f\n", total_production, g_zone_vector[orig].gravity_production[at], percentage_difference);

            }

//            fprintf(g_pFileODMatrix, "est,");

            //for (int dest = 0; dest < g_zone_vector.size(); ++dest)  // d
            //{
            //    fprintf(g_pFileODMatrix, "%f,", g_zone_vector[dest].gravity_est_attraction[at]);
            //}
            //fprintf(g_pFileODMatrix, "\n");
            //fprintf(g_pFileODMatrix, "target,");

            //for (int dest = 0; dest < g_zone_vector.size(); ++dest)  // d
            //{
            //    fprintf(g_pFileODMatrix, "%f,", g_zone_vector[dest].gravity_attraction[at]);
            //}

            //fprintf(g_pFileODMatrix, "\n");
            //fprintf(g_pFileODMatrix, "ratio_diff,");

            //for (int dest = 0; dest < g_zone_vector.size(); ++dest)  // d
            //{
            //    float percentage_difference = (g_zone_vector[dest].gravity_est_attraction[at] - g_zone_vector[dest].gravity_attraction[at]) / max(0.001, g_zone_vector[dest].gravity_attraction[at]);

            //    fprintf(g_pFileODMatrix, "%f,", percentage_difference);
            //}
            fprintf(g_pFileODMatrix, "\n");

            break; // only for one agent type
        }
        break; // only for one demand period
    }


    fclose(g_pFileODMatrix);
}


////////////////////set reference half of link capacity as link volume

//for (int l = 0; l < g_link_vector.size(); ++l)
//{
//    if(g_link_vector[l].lane_capacity >=800 && g_link_vector[l].lane_capacity < 2000) // remark: link_type = -1 is virtual connector
//    {
//        g_link_vector[l].obs_count = g_link_vector[l].lane_capacity * g_link_vector[l].number_of_lanes/2;
//    }
//}

for (int z = 0; z < g_zone_vector.size(); ++z)  // d
{
    g_zone_vector[z].obs_production = -1;
    g_zone_vector[z].obs_attraction = -1;  // invalidate the data, we will focus on link count data first 
}
}
void g_demand_file_generation(Assignment& assignment)
{
        g_trip_generation(assignment);
        g_writing_demand_files(assignment);
  
    
}