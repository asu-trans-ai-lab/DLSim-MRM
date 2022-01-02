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


void g_add_new_access_link(int internal_from_node_seq_no, int internal_to_node_seq_no, float link_distance_in_km, int agent_type_no, int zone_seq_no = -1)
{
    // create a link object
    CLink link;

    link.b_automated_generated_flag = true;
    link.from_node_seq_no = internal_from_node_seq_no;
    link.to_node_seq_no = internal_to_node_seq_no;
    link.link_seq_no = assignment.g_number_of_links;
    link.to_node_seq_no = internal_to_node_seq_no;

    link.link_type = 1000;  // access_link

    //only for outgoing connectors
    link.zone_seq_no_for_outgoing_connector = zone_seq_no;

    link.link_type_code = "access_link";
    //BPR
    link.traffic_flow_code = 0;

    link.spatial_capacity_in_vehicles = 99999;
    link.lane_capacity = 999999;
    link.link_spatial_capacity = 99999;
    link.link_distance_in_km = link_distance_in_km;
    link.free_speed = assignment.g_AgentTypeVector[agent_type_no].access_speed;

    for (int tau = 0; tau < assignment.g_number_of_demand_periods; ++tau)
    {
        //setup default values
        link.VDF_period[tau].period_capacity = 99999;
        // 60.0 for 60 min per hour
        link.fftt = link_distance_in_km / max(0.001, link.free_speed) * 60; 
        link.VDF_period[tau].FFTT = link.fftt;
        link.VDF_period[tau].penalty = 99;
        link.VDF_period[tau].alpha = 0;
        link.VDF_period[tau].beta = 0;
        link.VDF_period[tau].allowed_uses += assignment.g_AgentTypeVector[agent_type_no].agent_type;
        link.TDBaseTT[tau] = link.fftt;
        link.TDBaseCap[tau] = 99999;
        link.travel_time_per_period[tau] = link.fftt;
    }

    // add this link to the corresponding node as part of outgoing node/link
    g_node_vector[internal_from_node_seq_no].m_outgoing_link_seq_no_vector.push_back(link.link_seq_no);
    // add this link to the corresponding node as part of outgoing node/link
    g_node_vector[internal_to_node_seq_no].m_incoming_link_seq_no_vector.push_back(link.link_seq_no);
    // add this link to the corresponding node as part of outgoing node/link
    g_node_vector[internal_from_node_seq_no].m_to_node_seq_no_vector.push_back(link.to_node_seq_no);
    // add this link to the corresponding node as part of outgoing node/link
    g_node_vector[internal_from_node_seq_no].m_to_node_2_link_seq_no_map[link.to_node_seq_no] = link.link_seq_no;

    g_link_vector.push_back(link);

    assignment.g_number_of_links++;
}

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
                    double near_by_distance = g_calculate_p2p_distance_in_meter_from_latitude_longitude(g_node_vector[i].x, g_node_vector[i].y, g_node_vector[j].x, g_node_vector[j].y);

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

    FILE* g_pFileZone = nullptr;
    g_pFileZone = fopen("zone.csv", "w");

    if (g_pFileZone == NULL)
    {
        cout << "File zone.csv cannot be opened." << endl;
        g_program_stop();
    }


        fprintf(g_pFileZone, "node_id,zone_id,access_node_vector,cell_code,cell_id,access_distance,x_coord,y_coord,");
        fprintf(g_pFileZone, "geometry,");
        for (int at = 0; at < assignment.g_AgentTypeVector.size(); ++at)
        {
            fprintf(g_pFileZone, "%s_production,%s_attraction,", assignment.g_AgentTypeVector[at].agent_type.c_str(), assignment.g_AgentTypeVector[at].agent_type.c_str());
        }

        fprintf(g_pFileZone, "\n");


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
                __int64 cell_id = g_get_cell_ID(g_node_vector[i].x, g_node_vector[i].y, assignment.m_GridResolution);
                int zone_id;

                if (assignment.cell_id_mapping.find(cell_id) == assignment.cell_id_mapping.end())  // create a cell
                {
                    //create zone
                    assignment.cell_id_mapping[cell_id] = g_node_vector[i].node_id;
                    string cell_code = g_get_cell_code(g_node_vector[i].x, g_node_vector[i].y, assignment.m_GridResolution, left, top);


                    int x_i = floor(g_node_vector[i].x / assignment.m_GridResolution);
                    int y_i = floor(g_node_vector[i].y/ assignment.m_GridResolution);

                    double x_coord_left = x_i * assignment.m_GridResolution;
                    double y_coord_bottom = y_i * assignment.m_GridResolution;
                    double x_coord_right = x_coord_left + assignment.m_GridResolution;
                    double y_coord_top = y_coord_bottom + assignment.m_GridResolution;

                    fprintf(g_pFileZone, "%d,", 100000+assignment.cell_id_mapping.size()+1);
                    fprintf(g_pFileZone, "%d,", assignment.cell_id_mapping.size()+1);
                    // generate access nodes
                    std::vector <int> access_node_vector;
                    float max_distance = 0;

                    float zone_x = (g_node_vector[i].x * 2 + x_coord_left + x_coord_right) / 4;
                    float zone_y = (g_node_vector[i].y * 2 + y_coord_top + y_coord_bottom) / 4;
                    for (int j = 0; j < g_node_vector.size(); j++)
                    {

                        if (g_node_vector[j].is_activity_node >= 1)
                        {

                            __int64 cell_id_j = g_get_cell_ID(g_node_vector[j].x, g_node_vector[j].y, assignment.m_GridResolution);
                            if (cell_id == cell_id_j)
                            {
                                double distance = g_calculate_p2p_distance_in_meter_from_latitude_longitude(zone_x, zone_y, g_node_vector[j].x, g_node_vector[j].y);
                                if (distance > max_distance)
                                {
                                    max_distance = distance;
                                }

                                access_node_vector.push_back(g_node_vector[j].node_id);
                                fprintf(g_pFileZone, "%d;", g_node_vector[j].node_id);
                            }
                        }
                    }

                    fprintf(g_pFileZone, ",%s,%jd,", cell_code.c_str(), assignment.cell_id_mapping[cell_id]);
                    fprintf(g_pFileZone, "%f,", max_distance);


                    fprintf(g_pFileZone, "%f,%f,", (g_node_vector[i].x*2 + x_coord_left+ x_coord_right)/4, (g_node_vector[i].y * 2 + y_coord_top+ y_coord_bottom)/4);

                    fprintf(g_pFileZone, "\"LINESTRING (");

                    fprintf(g_pFileZone, "%f %f,", x_coord_left, y_coord_top);
                    fprintf(g_pFileZone, "%f %f,", x_coord_right, y_coord_top);
                    fprintf(g_pFileZone, "%f %f,", x_coord_right, y_coord_bottom);
                    fprintf(g_pFileZone, "%f %f,", x_coord_left, y_coord_bottom);
                    fprintf(g_pFileZone, "%f %f,", x_coord_left, y_coord_top);
                    fprintf(g_pFileZone, ")\",");

                    for (int at = 0; at < assignment.g_AgentTypeVector.size(); ++at)
                    {
                        fprintf(g_pFileZone, "0,0,");
                    }


                    fprintf(g_pFileZone, "\n");

                }




            }
        }

        dtalog.output() << "Step 1.4.3: creating " << assignment.cell_id_mapping.size() << " zones." << endl;
        fclose(g_pFileZone);
    
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
        ozone.cell_code = assignment.cell_id_2_cell_code_mapping[ozone.cell_id];
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
                        float distance_in_mile = g_calculate_p2p_distance_in_meter_from_latitude_longitude(g_zone_vector[orig].cell_x, g_zone_vector[orig].cell_y, g_zone_vector[dest].cell_x, g_zone_vector[dest].cell_y);
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
// 
//}

for (int z = 0; z < g_zone_vector.size(); ++z)  // d
{
    g_zone_vector[z].obs_production = -1;
    g_zone_vector[z].obs_attraction = -1;  // invalidate the data, we will focus on link count data first 
}
}
void g_demand_file_generation(Assignment& assignment)
{

//        g_trip_generation(assignment);
//        g_writing_demand_files(assignment);
  
}


void g_zone_to_access(Assignment& assignment)
{
    int debug_line_count = 0;
    if (assignment.assignment_mode == 21 && assignment.g_AgentTypeVector.size() > 0)  //zone2connector
    {
        int at = 0;

            if (assignment.g_AgentTypeVector[at].access_node_type.size() == 0)  // for each simple agent type
            {
                // find the closest zone id

                if (debug_line_count <= 20)
                {

                    dtalog.output() << " connector generation condition 1: agent type " << assignment.g_AgentTypeVector[at].agent_type.c_str() << " has access node type" << assignment.g_AgentTypeVector[at].access_node_type.size() << endl;
                    // zone without multimodal access
                    debug_line_count++;
                }

                for (int a_k = 0; a_k < g_node_vector.size(); a_k++)
                {
                    if (g_node_vector[a_k].is_activity_node == 100) //first loop for mode_specific activity node with zone id 
                    {

                        int zone_id = g_node_vector[a_k].zone_org_id;
                        int zone_seq_no = zone_seq_no = assignment.g_zoneid_to_zone_seq_no_mapping[zone_id];
                        float access_distance = g_node_vector[a_k].access_distance;

                        if (debug_line_count <= 20)
                        {

                            dtalog.output() << " connector generation generation condition 2: agent type no = " << at << " for node no. " << a_k << "as activity node with zone_id >=1" << endl;
                            // zone without multimodal access
                            debug_line_count++;
                        }


                        // stage 2:  // min_distance
                        double min_distance = 9999999;
                        int min_distance_node_seq_no = -1;

                        // stage 1:  // preferreed distance range
                        double min_distance_within_range = 9999999;
                        int min_distance_node_id_within_range = -1;

                        std::vector<int> access_node_seq_vector;
                        std::vector<float> access_node_distance_vector;

                        for (int i = 0; i < g_node_vector.size(); i++)
                        {
                            if (g_node_vector[i].is_activity_node == 2)  // is boundary
                            {

                                if (assignment.g_AgentTypeVector[at].access_node_type.find(g_node_vector[i].node_type) != string::npos)  // check allowed access code
                                {

                                    double zone_x = g_node_vector[a_k].x;
                                    double zone_y = g_node_vector[a_k].y;

                                    //test                                double near_by_distance_1 = g_calculate_p2p_distance_in_meter_from_latitude_longitude(-77.429293, 39.697895, -77.339847, 38.947676);

                                    double distance = g_calculate_p2p_distance_in_meter_from_latitude_longitude(zone_x, zone_y, g_node_vector[i].x, g_node_vector[i].y);
                                    // calculate the distance 

                                    if (distance < min_distance)
                                    {
                                        min_distance = distance;
                                        min_distance_node_seq_no = i;
                                    }

                                    if (distance <= access_distance*1.01)  // check the range 
                                    {
                                        min_distance_within_range = distance;
                                        min_distance_node_id_within_range = i;
                                        access_node_seq_vector.push_back(i);
                                        access_node_distance_vector.push_back(distance);
                                    }
                                }

                            }
                        }  // scan for all nodes


                        // check access node vector for each pair of zone and agent type
                        // 
                        if (access_node_seq_vector.size() > 0)  // preferred: access link within the range 
                        {


                            float distance_k_cut_off_value = 99999;

                            int acecss_link_k = 4;

                            if (access_node_distance_vector.size() > acecss_link_k)
                            {

                                std::vector<float> access_node_distance_vector_temp;
                                access_node_distance_vector_temp = access_node_distance_vector;
                                std::sort(access_node_distance_vector_temp.begin(), access_node_distance_vector_temp.end());

                                distance_k_cut_off_value = access_node_distance_vector_temp[max(0, assignment.g_AgentTypeVector[at].acecss_link_k - 1)];
                                //distance_k can be dynamically determined based on the density of stops and stations at different areas, e.g.CBM vs. rual area
                            }

                            for (int an = 0; an < access_node_seq_vector.size(); an++)
                            {
                                if (access_node_distance_vector[an] < distance_k_cut_off_value)  // within the shortest k ranage 
                                {
                                    if (g_node_vector[access_node_seq_vector[an]].is_boundary == 1)
                                        g_add_new_access_link(a_k, access_node_seq_vector[an], access_node_distance_vector[an], at, -1);

                                    if (g_node_vector[access_node_seq_vector[an]].is_boundary == -1)
                                        g_add_new_access_link(access_node_seq_vector[an], a_k, access_node_distance_vector[an], at, -1);

                                    assignment.g_AgentTypeVector[at].zone_id_cover_map[zone_id] = true;
                                }
                            }

                        }

                    }  // for each zone

                }

        }// for each agent type 

        g_OutputModelFiles(3);  // node
        g_program_exit();

    }
}