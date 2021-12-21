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

void g_load_scenario_data(Assignment& assignment)
{
    dtalog.output() << "Step 2.0: Reading scenario data..." << endl;

    dtalog.output() << "Step 1.8: Reading file section [demand_file_list] in setting.csv..." << endl;

    CCSVParser parser;
    parser.IsFirstLineHeader = false;

    int count = 0;

    if (parser.OpenCSVFile("scenario.csv", false))
    {
        while (parser.ReadRecord_Section())
        {
            if (parser.SectionName == "[capacity_scenario]")
            {
                int from_node_id = 0;
                if (!parser.GetValueByFieldName("from_node_id", from_node_id))
                {
                    dtalog.output() << "Error: from_node_id in file scenario.csv is not defined." << endl;
                    continue;
                }

                int to_node_id = 0;
                if (!parser.GetValueByFieldName("to_node_id", to_node_id))
                    continue;


                if (from_node_id == 22286 && to_node_id == 22290)
                { 
                    int idebug = 1;
                }

                if (assignment.g_node_id_to_seq_no_map.find(from_node_id) == assignment.g_node_id_to_seq_no_map.end())
                {
                    dtalog.output() << "Error: from_node_id " << from_node_id << " in file scenario.csv is not defined in node.csv." << endl;
                    //has not been defined
                    continue;
                }
                if (assignment.g_node_id_to_seq_no_map.find(to_node_id) == assignment.g_node_id_to_seq_no_map.end())
                {
                    dtalog.output() << "Error: to_node_id " << to_node_id << " in file scenario.csv is not defined in node.csv." << endl;
                    //has not been defined
                    continue;
                }

                // create a link object


                // map external node number to internal node seq no.
                int internal_from_node_seq_no = assignment.g_node_id_to_seq_no_map[from_node_id];
                int internal_to_node_seq_no = assignment.g_node_id_to_seq_no_map[to_node_id];

                int link_seq_no = 0;
                if (g_node_vector[internal_from_node_seq_no].m_to_node_2_link_seq_no_map.find(internal_to_node_seq_no) != g_node_vector[internal_from_node_seq_no].m_to_node_2_link_seq_no_map.end())
                {
                    link_seq_no = g_node_vector[internal_from_node_seq_no].m_to_node_2_link_seq_no_map[internal_to_node_seq_no];
                }
                else
                {
                    dtalog.output() << "Error: Link " << from_node_id << "->" << to_node_id << " in file scenario.csv is not defined in link.csv." << endl;
                    continue;
                }

                string time_period;
                if (!parser.GetValueByFieldName("time_period", time_period))
                {
                    dtalog.output() << "Error: Field time_window in file scenario.csv cannot be read." << endl;
                    g_program_stop();
                    break;
                }

                vector<float> global_minute_vector;

                //input_string includes the start and end time of a time period with hhmm format
                global_minute_vector = g_time_parser(time_period); //global_minute_vector incldue the starting and ending time
                if (global_minute_vector.size() == 2)
                {
                    if (global_minute_vector[0] < assignment.g_LoadingStartTimeInMin)
                        global_minute_vector[0] = assignment.g_LoadingStartTimeInMin;

                    if (global_minute_vector[0] > assignment.g_LoadingEndTimeInMin)
                        global_minute_vector[0] = assignment.g_LoadingEndTimeInMin;

                    if (global_minute_vector[1] < assignment.g_LoadingStartTimeInMin)
                        global_minute_vector[1] = assignment.g_LoadingStartTimeInMin;

                    if (global_minute_vector[1] > assignment.g_LoadingEndTimeInMin)
                        global_minute_vector[1] = assignment.g_LoadingEndTimeInMin;

                    if (global_minute_vector[1] < global_minute_vector[0])
                        global_minute_vector[1] = global_minute_vector[0];

                }
                else
                    continue;

                // capacity in the space time arcs
                float capacity = 1;
                parser.GetValueByFieldName("capacity", capacity);


                string demand_period;

                parser.GetValueByFieldName("demand_period", demand_period);


                unsigned int RandomSeed = 101;
                float residual;
                float random_ratio = 0;

                for (int s = global_minute_vector[0] *60; s <= global_minute_vector[1]*60; s++)
                {
                    int t_simu_second = (s- assignment.g_LoadingStartTimeInMin * 60);

                    if (t_simu_second < 0)
                        t_simu_second = 0;
                    double per_sec_capacity = capacity / 3600;

                    if (capacity < 1)
                    {
                        g_link_vector[link_seq_no].m_link_pedefined_capacity_map[t_simu_second] = 0;
                    }
                    else
                    {
                        residual = per_sec_capacity - (int)(per_sec_capacity);
                        RandomSeed = (LCG_a * RandomSeed + LCG_c) % LCG_M;
                        random_ratio = float(RandomSeed) / LCG_M;

                        if (random_ratio < residual)
                            g_link_vector[link_seq_no].m_link_pedefined_capacity_map[t_simu_second] = 1;
                        else
                            g_link_vector[link_seq_no].m_link_pedefined_capacity_map[t_simu_second] = 0;
                    }


                }


                int tau = 0;
                if (assignment.demand_period_to_seqno_mapping.find(demand_period) != assignment.demand_period_to_seqno_mapping.end())
                {
                    tau = assignment.demand_period_to_seqno_mapping[demand_period];
                }
                //else
                //{
                //    dtalog.output() << "demand_period= " << demand_period.c_str() << " in scenario.csv has not been defined in setting.csv." << endl;
                //    g_program_stop();
                //}


                //
                if (capacity < g_link_vector[link_seq_no].lane_capacity * g_link_vector[link_seq_no].number_of_lanes * 0.8)
                {
                    g_link_vector[link_seq_no].capacity_reduction_map[tau] = 1;
                }

                char lr_price_field_name[50];

                for (int at = 0; at < assignment.g_AgentTypeVector.size(); at++)
                {
                    double LR_ini_price = 0;
                    double LR_RT_price = 0;

                    sprintf(lr_price_field_name, "lr_price_%s", assignment.g_AgentTypeVector[at].agent_type.c_str());
                    parser.GetValueByFieldName(lr_price_field_name, LR_ini_price,true,false);
                    g_link_vector[link_seq_no].VDF_period[tau].LR_price[at] = LR_ini_price;

                    if(assignment.g_AgentTypeVector[at].real_time_information>=1)
                    {
                    sprintf(lr_price_field_name, "lr_rt_price_%s", assignment.g_AgentTypeVector[at].agent_type.c_str());
                    parser.GetValueByFieldName(lr_price_field_name, LR_RT_price, true, false);
                    g_link_vector[link_seq_no].VDF_period[tau].LR_RT_price[at] = LR_RT_price;

                        if (fabs(LR_RT_price) > 0.001)
                        {
                            dtalog.output() << "link " << from_node_id << "->" << to_node_id << " has a lr RT price of " << g_link_vector[link_seq_no].VDF_period[tau].LR_RT_price[at] << " for agent type "
                                << assignment.g_AgentTypeVector[at].agent_type.c_str() << " at demand period " << demand_period.c_str() << endl;
                        }

                    }

                    if (fabs(LR_ini_price) > 0.001)
                    {
                        dtalog.output() << "link " << from_node_id << "->" << to_node_id << " has a lr price of " << g_link_vector[link_seq_no].VDF_period[tau].LR_price[at] << " for agent type "
                            << assignment.g_AgentTypeVector[at].agent_type.c_str() << " at demand period " << demand_period.c_str() << endl;
                    }

                }

                count++;

               dtalog.output() << "reading " << count << " capacity reduction scenario.. " << endl;
            }
        }

        parser.CloseCSVFile();
    }
    // we now know the number of links

}
