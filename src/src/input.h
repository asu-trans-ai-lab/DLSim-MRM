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

void g_ReadDemandFileBasedOnDemandFileList(Assignment& assignment)
{
    //	fprintf(g_pFileOutputLog, "number of zones =,%lu\n", g_zone_vector.size());

    assignment.InitializeDemandMatrix(g_zone_vector.size(), assignment.g_AgentTypeVector.size(), assignment.g_DemandPeriodVector.size());

    float total_demand_in_demand_file = 0;

    CCSVParser parser;
    dtalog.output() << endl;
    dtalog.output() << "Step 1.8: Reading file section [demand_file_list] in setting.csv..." << endl;
    parser.IsFirstLineHeader = false;

    if (parser.OpenCSVFile("settings.csv", false))
    {
        while (parser.ReadRecord_Section())
        {
            if (parser.SectionName == "[demand_file_list]")
            {
                int file_sequence_no = 1;

                string format_type = "null";

                int demand_format_flag = 0;

                if (!parser.GetValueByFieldName("file_sequence_no", file_sequence_no))
                    break;

                // skip negative sequence no
                if (file_sequence_no <= -1)
                    continue;

                double loading_scale_factor = 1.0;
                string file_name, demand_period, agent_type;
                parser.GetValueByFieldName("file_name", file_name);
                parser.GetValueByFieldName("demand_period", demand_period);
                parser.GetValueByFieldName("format_type", format_type);
                parser.GetValueByFieldName("loading_scale_factor", loading_scale_factor, false);



                if (format_type.find("null") != string::npos)  // skip negative sequence no
                {
                    dtalog.output() << "Please provide format_type in section [demand_file_list.]" << endl;
                    g_ProgramStop();
                }

                parser.GetValueByFieldName("agent_type", agent_type);

                int agent_type_no = 0;
                int demand_period_no = 0;

                if (assignment.demand_period_to_seqno_mapping.find(demand_period) != assignment.demand_period_to_seqno_mapping.end())
                    demand_period_no = assignment.demand_period_to_seqno_mapping[demand_period];
                else
                {
                    dtalog.output() << "Error: demand period in section [demand_file_list]" << demand_period << " cannot be found." << endl;
                    g_ProgramStop();
                }

                bool b_multi_agent_list = false;

                if (agent_type == "multi_agent_list")
                    b_multi_agent_list = true;
                else
                {
                    if (assignment.agent_type_2_seqno_mapping.find(agent_type) != assignment.agent_type_2_seqno_mapping.end())
                        agent_type_no = assignment.agent_type_2_seqno_mapping[agent_type];
                    else
                    {
                        dtalog.output() << "Error: agent_type in agent_type " << agent_type << " cannot be found." << endl;
                        g_ProgramStop();
                    }
                }

                if (demand_period_no > _MAX_TIMEPERIODS)
                {
                    dtalog.output() << "demand_period_no should be less than settings in demand_period section. Please change the parameter settings in the source code." << endl;
                    g_ProgramStop();
                }

                if (format_type.find("column") != string::npos)  // or muliti-column
                {
                    bool bFileReady = false;
                    int error_count = 0;
                    int critical_OD_count = 0;
                    double critical_OD_volume = 0;

                    // read the file formaly after the test.
                    FILE* st;
                    fopen_ss(&st, file_name.c_str(), "r");
                    if (st)
                    {
                        bFileReady = true;
                        int line_no = 0;

                        dtalog.output() << endl << "VIODP  o,d,volume,geometry" << endl;

                        while (true)
                        {
                            int origin_zone = (int)(g_read_float(st));
                            int destination_zone = (int)g_read_float(st);
                            float demand_value = g_read_float(st);

                            if (origin_zone <= -1)
                            {
                                if (line_no == 1 && !feof(st))  // read only one line, but has not reached the end of the line
                                {
                                    dtalog.output() << endl << "Error: Only one line has been read from file. Are there multiple columns of demand type in file " << file_name << " per line?" << endl;
                                    g_ProgramStop();
                                }
                                break;
                            }

                            if (assignment.g_zoneid_to_zone_seq_no_mapping.find(origin_zone) == assignment.g_zoneid_to_zone_seq_no_mapping.end())
                            {
                                if (error_count < 10)
                                    dtalog.output() << endl << "Warning: origin zone " << origin_zone << "  has not been defined in node.csv" << endl;

                                error_count++;
                                // origin zone has not been defined, skipped.
                                continue;
                            }

                            if (assignment.g_zoneid_to_zone_seq_no_mapping.find(destination_zone) == assignment.g_zoneid_to_zone_seq_no_mapping.end())
                            {
                                if (error_count < 10)
                                    dtalog.output() << endl << "Warning: destination zone " << destination_zone << "  has not been defined in node.csv" << endl;

                                error_count++;
                                // destination zone has not been defined, skipped.
                                continue;
                            }

                            int from_zone_seq_no = 0;
                            int to_zone_seq_no = 0;
                            from_zone_seq_no = assignment.g_zoneid_to_zone_seq_no_mapping[origin_zone];
                            to_zone_seq_no = assignment.g_zoneid_to_zone_seq_no_mapping[destination_zone];

                            // encounter return
                            if (demand_value < -99)
                                break;

                            demand_value *= loading_scale_factor;
                            if (demand_value >= 5)
                            {
                                critical_OD_volume += demand_value;
                                critical_OD_count += 1;
                                //dtalog.output() << origin_zone << "," << destination_zone << "," << demand_value << "," << "\"LINESTRING( " <<
                                //    assignment.zone_id_X_mapping[origin_zone] << " " << assignment.zone_id_Y_mapping[origin_zone] << "," <<
                                //    assignment.zone_id_X_mapping[destination_zone] << " " << assignment.zone_id_Y_mapping[destination_zone] << ")\" " << endl;

                            }

                            assignment.total_demand[agent_type_no][demand_period_no] += demand_value;
                            assignment.g_column_pool[from_zone_seq_no][to_zone_seq_no][agent_type_no][demand_period_no].od_volume += demand_value;
                            assignment.total_demand_volume += demand_value;
                            assignment.g_origin_demand_array[from_zone_seq_no] += demand_value;

                            // we generate vehicles here for each OD data line
                            if (line_no <= 5)  // read only one line, but has not reached the end of the line
                                dtalog.output() << "o_zone_id:" << origin_zone << ", d_zone_id: " << destination_zone << ", value = " << demand_value << endl;

                            line_no++;
                        }  // scan lines

                        fclose(st);

                        dtalog.output() << "total demand volume is " << assignment.total_demand_volume << endl;
                        dtalog.output() << "crtical demand volume has " << critical_OD_count << " OD pairs in size," << critical_OD_volume << ", " << ", account for " << critical_OD_volume / max(0.1, assignment.total_demand_volume) * 100 << "%%" << endl;

                        dtalog.output() << "crtical OD zones volume has " << critical_OD_count << " OD pairs in size," << critical_OD_volume << ", " << ", account for " << critical_OD_volume / max(0.1, assignment.total_demand_volume) * 100 << "%%" << endl;


                        std::map<int, float>::iterator it;
                        int count_zone_demand = 0;
                        for (it = assignment.g_origin_demand_array.begin(); it != assignment.g_origin_demand_array.end(); ++it)
                        {
                            //if (it->second > 5)
                            //{
                            //    dtalog.output() << "o_zone " << it->first << ", d_zone=," << it->second << endl;
                            //    count_zone_demand++;
                            //}
                        }
                        dtalog.output() << "There are  " << count_zone_demand << " zones with positive demand" << endl;

                    }
                    else
                    {
                        // open file
                        dtalog.output() << "Error: File " << file_name << " cannot be opened.\n It might be currently used and locked by EXCEL." << endl;
                        g_ProgramStop();
                    }
                }
                else if (format_type.compare("path") == 0)
                {

                    int path_counts = 0;
                    float sum_of_path_volume = 0;
                    CCSVParser parser;
                    if (parser.OpenCSVFile(file_name, false))
                    {
                        int total_path_in_demand_file = 0;
                        // read agent file line by line,

                        int agent_id, o_zone_id, d_zone_id;
                        string agent_type, demand_period;

                        std::vector <int> node_sequence;

                        while (parser.ReadRecord())
                        {
                            total_path_in_demand_file++;
                            if (total_path_in_demand_file % 1000 == 0)
                                dtalog.output() << "total_path_in_demand_file is " << total_path_in_demand_file << endl;

                            parser.GetValueByFieldName("agent_id", agent_id);
                            parser.GetValueByFieldName("o_zone_id", o_zone_id);
                            parser.GetValueByFieldName("d_zone_id", d_zone_id);

                            CAgentPath agent_path_element;

                            agent_path_element.path_id = 0;
                            parser.GetValueByFieldName("path_id", agent_path_element.path_id, false);


                            int from_zone_seq_no = 0;
                            int to_zone_seq_no = 0;
                            from_zone_seq_no = assignment.g_zoneid_to_zone_seq_no_mapping[o_zone_id];
                            to_zone_seq_no = assignment.g_zoneid_to_zone_seq_no_mapping[d_zone_id];

                            if (format_type.compare("path") == 0)
                            {
                                double volume = 0;
                                parser.GetValueByFieldName("volume", volume);
                                volume *= loading_scale_factor;
                                agent_path_element.volume = volume;
                                path_counts++;
                                sum_of_path_volume += agent_path_element.volume;

                                assignment.total_demand[agent_type_no][demand_period_no] += agent_path_element.volume;
                                assignment.g_column_pool[from_zone_seq_no][to_zone_seq_no][agent_type_no][demand_period_no].od_volume += agent_path_element.volume;
                                assignment.total_demand_volume += agent_path_element.volume;
                                assignment.g_origin_demand_array[from_zone_seq_no] += agent_path_element.volume;
                            }



                            //apply for both agent csv and routing policy
                            assignment.g_column_pool[from_zone_seq_no][to_zone_seq_no][agent_type_no][demand_period_no].bfixed_route = true;

                            bool bValid = true;

                            string path_node_sequence;
                            parser.GetValueByFieldName("node_sequence", path_node_sequence);

                            std::vector<int> node_id_sequence;

                            g_ParserIntSequence(path_node_sequence, node_id_sequence);

                            std::vector<int> node_no_sequence;
                            std::vector<int> link_no_sequence;

                            int node_sum = 0;
                            for (int i = 0; i < node_id_sequence.size(); ++i)
                            {
                                if (assignment.g_node_id_to_seq_no_map.find(node_id_sequence[i]) == assignment.g_node_id_to_seq_no_map.end())
                                {
                                    bValid = false;
                                    //has not been defined
                                    continue;
                                    // warning
                                }

                                int internal_node_seq_no = assignment.g_node_id_to_seq_no_map[node_id_sequence[i]];  // map external node number to internal node seq no.
                                node_no_sequence.push_back(internal_node_seq_no);

                                node_sum += internal_node_seq_no;
                                if (i >= 1)
                                {
                                    // check if a link exists
                                    int link_seq_no = -1;
                                    // map external node number to internal node seq no.
                                    int prev_node_seq_no = assignment.g_node_id_to_seq_no_map[node_id_sequence[i - 1]];
                                    int current_node_no = node_no_sequence[i];

                                    if (g_node_vector[prev_node_seq_no].m_to_node_2_link_seq_no_map.find(current_node_no) != g_node_vector[prev_node_seq_no].m_to_node_2_link_seq_no_map.end())
                                    {
                                        link_seq_no = g_node_vector[prev_node_seq_no].m_to_node_2_link_seq_no_map[node_no_sequence[i]];
                                        link_no_sequence.push_back(link_seq_no);
                                    }
                                    else
                                        bValid = false;
                                }
                            }

                            if (bValid)
                            {
                                agent_path_element.node_sum = node_sum; // pointer to the node sum based path node sequence;
                                agent_path_element.path_link_sequence = link_no_sequence;

                                CColumnVector* pColumnVector = &(assignment.g_column_pool[from_zone_seq_no][to_zone_seq_no][agent_type_no][demand_period_no]);

                                // we cannot find a path with the same node sum, so we need to add this path into the map,
                                if (pColumnVector->path_node_sequence_map.find(node_sum) == pColumnVector->path_node_sequence_map.end())
                                {
                                    // add this unique path
                                    int path_count = pColumnVector->path_node_sequence_map.size();
                                    pColumnVector->path_node_sequence_map[node_sum].path_seq_no = path_count;
                                    pColumnVector->path_node_sequence_map[node_sum].path_volume = 0;
                                    //assignment.g_column_pool[m_origin_zone_seq_no][destination_zone_seq_no][agent_type][tau].time = m_label_time_array[i];
                                    //assignment.g_column_pool[m_origin_zone_seq_no][destination_zone_seq_no][agent_type][tau].path_node_sequence_map[node_sum].path_distance = m_label_distance_array[i];
                                    pColumnVector->path_node_sequence_map[node_sum].path_toll = 0;

                                    pColumnVector->path_node_sequence_map[node_sum].AllocateVector(node_no_sequence, link_no_sequence, false);
                                }

                                pColumnVector->path_node_sequence_map[node_sum].path_volume += agent_path_element.volume;
                            }
                        }
                        dtalog.output() << "total_demand_volume loaded from path file is " << sum_of_path_volume << " with " << path_counts << "paths." << endl;

                    }
                    else
                    {
                        //open file
                        dtalog.output() << "Error: File " << file_name << " cannot be opened.\n It might be currently used and locked by EXCEL." << endl;
                        g_ProgramStop();
                    }
                }
                else if (format_type.compare("activity_plan") == 0)
                {
                    /////////////////////
                    int path_counts = 0;
                    float sum_of_path_volume = 0;
                    CCSVParser parser;
                    if (parser.OpenCSVFile(file_name, false))
                    {
                        int total_path_in_demand_file = 0;
                        // read agent file line by line,

                        int agent_id, o_zone_id, d_zone_id;
                        string agent_type, demand_period;

                        std::vector <int> node_sequence;

                        while (parser.ReadRecord())
                        {
                            total_path_in_demand_file++;
                            if (total_path_in_demand_file % 1000 == 0)
                                dtalog.output() << "total_path_in_demand_file is " << total_path_in_demand_file << endl;

                            parser.GetValueByFieldName("agent_id", agent_id);
                            parser.GetValueByFieldName("o_zone_id", o_zone_id);
                            parser.GetValueByFieldName("d_zone_id", d_zone_id);

                            CAgentPath agent_path_element;


                            int from_zone_seq_no = 0;
                            int to_zone_seq_no = 0;  

                            // add protection later for activity and path data types
                            from_zone_seq_no = assignment.g_zoneid_to_zone_seq_no_mapping[o_zone_id];
                            to_zone_seq_no = assignment.g_zoneid_to_zone_seq_no_mapping[d_zone_id];

                            double volume = 0;
                            parser.GetValueByFieldName("volume", volume);
                            volume *= loading_scale_factor;
                            agent_path_element.volume = volume;
                            path_counts++;
                            sum_of_path_volume += agent_path_element.volume;

                            assignment.total_demand[agent_type_no][demand_period_no] += agent_path_element.volume;
                            assignment.g_column_pool[from_zone_seq_no][to_zone_seq_no][agent_type_no][demand_period_no].od_volume += agent_path_element.volume;
                            assignment.total_demand_volume += agent_path_element.volume;
                            assignment.g_origin_demand_array[from_zone_seq_no] += agent_path_element.volume;

                            //apply for both agent csv and routing policy

                            bool bValid = true;

                            string zone_sequence;
                            parser.GetValueByFieldName("zone_sequence", zone_sequence);

                            string agent_type_sequence;
                            parser.GetValueByFieldName("activity_agent_type", agent_type_sequence);
                            
                            

                            std::vector<int> zone_id_sequence;
                            std::vector<string> agent_type_string_sequence;

                            int activity_zone_size = g_ParserIntSequence(zone_sequence, zone_id_sequence);
                            int agent_type_size = g_ParserStringSequence(agent_type_sequence, agent_type_string_sequence);

                            //if (agent_type_size != activity_zone_size + 1)
                            //{
                            //    dtalog.output() << "Error: agent_type_size != activity_zone_size + 1" << endl;
                            //    g_ProgramStop();
                            //}


                            assignment.agent_type_2_seqno_mapping[agent_type_sequence];

                            std::vector<int> zone_no_sequence;

                            for (int i = 0; i < zone_id_sequence.size(); ++i)
                            {
                                if (assignment.g_zoneid_to_zone_seq_no_mapping.find(zone_id_sequence[i]) == assignment.g_zoneid_to_zone_seq_no_mapping.end())
                                {
                                    bValid = false;
                                    continue;
                                }

                                //if (assignment.agent_type_2_seqno_mapping.find(agent_type_string_sequence[i]) == assignment.agent_type_2_seqno_mapping.end())
                                //{
                                //    bValid = false;
                                //    continue;
                                //}



                                int internal_zone_seq_no = assignment.g_zoneid_to_zone_seq_no_mapping[zone_id_sequence[i]];  // map external node number to internal node seq no.
                                assignment.zone_seq_no_2_activity_mapping[internal_zone_seq_no] = 1; // notify the column pool to prepare the paths for this activity zones, even no direct OD volume 
                                zone_no_sequence.push_back(internal_zone_seq_no);
                            }

                            if (bValid)
                            {

                                CColumnVector* pColumnVector = &(assignment.g_column_pool[from_zone_seq_no][to_zone_seq_no][agent_type_no][demand_period_no]);
                                pColumnVector->activity_zone_no_vector.push_back(from_zone_seq_no);  // first origin
                                // to be protected later by checking
                                pColumnVector->activity_agent_type_no = assignment.agent_type_2_seqno_mapping[agent_type_sequence];

                                for(int i = 0; i < zone_no_sequence.size(); ++i)
                                {
                                pColumnVector->activity_zone_no_vector.push_back(zone_no_sequence[i]);  // for all intermedate activity zones
                                }

                                pColumnVector->activity_zone_no_vector.push_back(to_zone_seq_no);  // last destination 
                                // in total, we have # of activity zones + 2 as zone no sequence vector

                            }

                           }//end of file
                    } //open file
                    else
                    {
                        //open file
                        dtalog.output() << "Error: File " << file_name << " cannot be opened.\n It might be currently used and locked by EXCEL." << endl;
                        g_ProgramStop();
                    }

                }// activity
                else if (format_type.compare("matrix") == 0)
                {
                    bool bFileReady = false;
                    int error_count = 0;
                    int critical_OD_count = 0;
                    double critical_OD_volume = 0;

                    vector<int> LineIntegerVector;

                    CCSVParser parser;
                    parser.IsFirstLineHeader = false;
                    if (parser.OpenCSVFile(file_name, true))
                    {
                        int control_type_code;
                        int i = 0;
                        if (parser.ReadRecord())
                        {
                            parser.ConvertLineStringValueToIntegers();
                            LineIntegerVector = parser.LineIntegerVector;
                        }
                        parser.CloseCSVFile();
                    }

                    int number_of_zones = LineIntegerVector.size();


                    bFileReady = false;
                    int i;

                    FILE* st;
                    fopen_s(&st, file_name.c_str(), "r");
                    if (st != NULL)
                    {
                        // read the first line
                        g_read_a_line(st);

                        cout << "number of zones to be read = " << number_of_zones << endl;

                        //test if a zone has been defined. 
                        for (int destination_zone_index = 0; destination_zone_index < number_of_zones; destination_zone_index++)
                        {
                            int zone = LineIntegerVector[destination_zone_index];
                            if (assignment.g_zoneid_to_zone_seq_no_mapping.find(zone) == assignment.g_zoneid_to_zone_seq_no_mapping.end())
                            {
                                if (error_count < 10)
                                    dtalog.output() << endl << "Warning: destination zone " << zone << "  has not been defined in node.csv" << endl;

                                error_count++;
                                // destination zone has not been defined, skipped.
                                continue;
                            }

                        }


                        int line_no = 0;
                        for (int origin_zone_index = 0; origin_zone_index < number_of_zones; origin_zone_index++)
                        {
                            int origin_zone = (int)(g_read_float(st)); // read the origin zone number

                            if (assignment.g_zoneid_to_zone_seq_no_mapping.find(origin_zone) == assignment.g_zoneid_to_zone_seq_no_mapping.end())
                            {
                                if (error_count < 10)
                                    dtalog.output() << endl << "Warning: destination zone " << origin_zone << "  has not been defined in node.csv" << endl;

                                error_count++;
                                // destination zone has not been defined, skipped.
                                continue;
                            }

                            cout << "Reading file no." << file_sequence_no << " " << file_name << " at zone " << origin_zone << " ... " << endl;

                            for (int destination_zone_index = 0; destination_zone_index < number_of_zones; destination_zone_index++)
                            {
                                int destination_zone = LineIntegerVector[destination_zone_index];

                                float demand_value = g_read_float(st);

                                int from_zone_seq_no = 0;
                                int to_zone_seq_no = 0;
                                from_zone_seq_no = assignment.g_zoneid_to_zone_seq_no_mapping[origin_zone];
                                to_zone_seq_no = assignment.g_zoneid_to_zone_seq_no_mapping[destination_zone];

                                // encounter return
                                if (demand_value < -99)
                                    break;

                                demand_value *= loading_scale_factor;
                                if (demand_value >= 1)
                                {
                                    critical_OD_volume += demand_value;
                                    critical_OD_count += 1;
                                    //dtalog.output() << origin_zone << "," << destination_zone << "," << demand_value << "," << "\"LINESTRING( " <<
                                    //    assignment.zone_id_X_mapping[origin_zone] << " " << assignment.zone_id_Y_mapping[origin_zone] << "," <<
                                    //    assignment.zone_id_X_mapping[destination_zone] << " " << assignment.zone_id_Y_mapping[destination_zone] << ")\" " << endl;

                                }

                                assignment.total_demand[agent_type_no][demand_period_no] += demand_value;
                                assignment.g_column_pool[from_zone_seq_no][to_zone_seq_no][agent_type_no][demand_period_no].od_volume += demand_value;
                                assignment.total_demand_volume += demand_value;
                                assignment.g_origin_demand_array[from_zone_seq_no] += demand_value;

                                // we generate vehicles here for each OD data line
                                if (line_no <= 5)  // read only one line, but has not reached the end of the line
                                    dtalog.output() << "o_zone_id:" << origin_zone << ", d_zone_id: " << destination_zone << ", value = " << demand_value << endl;

                                line_no++;
                            }  // scan lines

                        }

                        fclose(st);

                        dtalog.output() << "total demand volume is " << assignment.total_demand_volume << endl;
                        dtalog.output() << "crtical demand volume has " << critical_OD_count << " OD pairs in size," << critical_OD_volume << ", " << ", account for " << critical_OD_volume / max(0.1, assignment.total_demand_volume) * 100 << "%%" << endl;

                        dtalog.output() << "crtical OD zones volume has " << critical_OD_count << " OD pairs in size," << critical_OD_volume << ", " << ", account for " << critical_OD_volume / max(0.1, assignment.total_demand_volume) * 100 << "%%" << endl;


                        std::map<int, float>::iterator it;
                        int count_zone_demand = 0;
                        for (it = assignment.g_origin_demand_array.begin(); it != assignment.g_origin_demand_array.end(); ++it)
                        {
                            if (it->second > 0.001)
                            {
                                dtalog.output() << "o_zone " << it->first << ", demand=," << it->second << endl;
                                count_zone_demand++;
                            }
                        }
                        dtalog.output() << "There are  " << count_zone_demand << " zones with positive demand" << endl;


                    } //end reading file
                    else
                    {
                        // open file
                        dtalog.output() << "Error: File " << file_name << " cannot be opened.\n It might be currently used and locked by EXCEL." << endl;
                        g_ProgramStop();
                    }

                }
                else {
                    dtalog.output() << "Error: format_type = " << format_type << " is not supported. Currently STALite supports format such as column, matrix, activity_plan, path." << endl;
                    g_ProgramStop();
                }
            }
        }
    }
}

void g_ReadOutputFileConfiguration(Assignment& assignment)
{
    dtalog.output() << "Step 1.9: Reading file section [output_file_configuration] in setting.csv..." << endl;

    cout << "Step 1.8: Reading file section [output_file_configuration] in setting.csv..." << endl;

    CCSVParser parser;
    parser.IsFirstLineHeader = false;
    if (parser.OpenCSVFile("settings.csv", false))
    {
        while (parser.ReadRecord_Section())
        {
            if (parser.SectionName == "[output_file_configuration]")
            {
                parser.GetValueByFieldName("path_output", assignment.path_output, false, false);
                parser.GetValueByFieldName("major_path_volume_threshold", assignment.major_path_volume_threshold, false, false);
                parser.GetValueByFieldName("trajectory_output", assignment.trajectory_output, false, false);
                parser.GetValueByFieldName("trajectory_sampling_rate", assignment.trajectory_sampling_rate, false, false);
                parser.GetValueByFieldName("trajectory_diversion_only", assignment.trajectory_diversion_only, false, false);
                parser.GetValueByFieldName("dynamic_link_performance_sampling_interval_in_min", assignment.dynamic_link_performance_sampling_interval_in_min, false, false);
                parser.GetValueByFieldName("dynamic_link_performance_sampling_interval_hd_in_min", assignment.dynamic_link_performance_sampling_interval_hd_in_min, false, false);

                dtalog.output() << "dynamic_link_performance_sampling_interval_in_min= " << assignment.dynamic_link_performance_sampling_interval_in_min << " min" << endl;
                dtalog.output() << "dynamic_link_performance_sampling_interval_hd_in_min= " << assignment.dynamic_link_performance_sampling_interval_hd_in_min << " min" << endl;

            }
        }

        parser.CloseCSVFile();
    }
}

void g_add_new_virtual_connector_link(int internal_from_node_seq_no, int internal_to_node_seq_no, string agent_type_str, int zone_seq_no = -1 )
{
    // create a link object
    CLink link;

    link.from_node_seq_no = internal_from_node_seq_no;
    link.to_node_seq_no = internal_to_node_seq_no;
    link.link_seq_no = assignment.g_number_of_links;
    link.to_node_seq_no = internal_to_node_seq_no;
    //virtual connector
    link.link_type = -1;
    

    //only for outgoing connectors
    link.zone_seq_no_for_outgoing_connector = zone_seq_no;

    //BPR
    link.traffic_flow_code = 0;

    link.spatial_capacity_in_vehicles = 99999;
    link.lane_capacity = 999999;
    link.link_spatial_capacity = 99999;
    link.length = 0.00001;

    for (int tau = 0; tau < assignment.g_number_of_demand_periods; ++tau)
    {
        //setup default values
        link.VDF_period[tau].capacity = 99999;
        // 60.0 for 60 min per hour
        link.VDF_period[tau].FFTT = 0;
        link.VDF_period[tau].alpha = 0;
        link.VDF_period[tau].beta = 0;
        link.VDF_period[tau].allowed_uses = agent_type_str;

        link.TDBaseTT[tau] = 0;
        link.TDBaseCap[tau] = 99999;
        link.travel_time_per_period[tau] = 0;

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

void g_add_new_access_link(int internal_from_node_seq_no, int internal_to_node_seq_no, float length, int agent_type_no)
{
    // create a link object
    CLink link;

    link.b_automated_generated_flag = true;
    link.from_node_seq_no = internal_from_node_seq_no;
    link.to_node_seq_no = internal_to_node_seq_no;
    link.link_seq_no = assignment.g_number_of_links;
    link.to_node_seq_no = internal_to_node_seq_no;

    link.link_type = 0;
    link.link_type_code = "multi-modal access_link";
    //BPR
    link.traffic_flow_code = 0;

    link.spatial_capacity_in_vehicles = 99999;
    link.lane_capacity = 999999;
    link.link_spatial_capacity = 99999;
    link.length = length;
    link.free_speed = assignment.g_AgentTypeVector[agent_type_no].access_speed;

    for (int tau = 0; tau < assignment.g_number_of_demand_periods; ++tau)
    {
        //setup default values
        link.VDF_period[tau].capacity = 99999;
        // 60.0 for 60 min per hour
        link.VDF_period[tau].FFTT = 0;
        link.VDF_period[tau].alpha = 0;
        link.VDF_period[tau].beta = 0;
        link.VDF_period[tau].allowed_uses += assignment.g_AgentTypeVector[agent_type_no].agent_type;
        link.TDBaseTT[tau] = 0;
        link.TDBaseCap[tau] = 99999;
        link.travel_time_per_period[tau] = 0;



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

double g_CalculateP2PDistanceInLonglatFromLatitudeLongitude(double p1x, double p1y, double p2x, double p2y)
{
    double distance = sqrt((p2y - p1y) * (p2y - p1y) + (p2x - p1x) * (p2x - p1x));
    return distance;
}


double g_CalculateP2PDistanceInMeterFromLatitudeLongitude(double p1x, double p1y, double p2x, double p2y)
{
    double PI = 3.1415926;
    double Equatorial_Radius = 3963.19059 * 1609; // unit: mile-> meter
    double toradians = 3.1415926 / 180.0;
    double todeg = 180.0 / PI;

    double p2lat = p2x * toradians;
    double p2lng = p2y * toradians;

    double p1lat = p1x * toradians;
    double p1lng = p1y * toradians;

    double distance = acos(sin(p1lat) * sin(p2lat) + cos(p1lat) * cos(p2lat) * cos(p2lng - p1lng)) * Equatorial_Radius;  // unit: mile
    return distance;
}

double g_CheckActivityNodes(Assignment& assignment)
{

    int activity_node_count = 0;
    for (int i = 0; i < g_node_vector.size(); i++)
    {

        if (g_node_vector[i].is_activity_node >= 1)
        {
            activity_node_count++;
        }
    }


    if (activity_node_count <= 1)
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

        //if (activity_node_count <= 1)
        //{
        //    activity_node_count = 0;
        //    sampling_rate = 2;

        //    for (int i = 0; i < g_node_vector.size(); i++)
        //    {

        //        if (i % sampling_rate == 0)
        //        {
        //            g_node_vector[i].is_activity_node = 10;//random generation
        //            activity_node_count++;
        //        }
        //    }
        //     still no activity nodes, define all nodes as activity nodes
        //    if (activity_node_count <= 1)
        //    {
        //        activity_node_count = 0;

        //        for (int i = 0; i < g_node_vector.size(); i++)
        //        {

        //            g_node_vector[i].is_activity_node = 10;//random generation
        //            activity_node_count++;
        //        }
        //    }
        //}


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
                    double near_by_distance = g_CalculateP2PDistanceInLonglatFromLatitudeLongitude(g_node_vector[i].x, g_node_vector[i].y, g_node_vector[j].x, g_node_vector[j].y);

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


void g_ReadInputData(Assignment& assignment)
{
    assignment.g_LoadingStartTimeInMin = 99999;
    assignment.g_LoadingEndTimeInMin = 0;

    //step 0:read demand period file
    CCSVParser parser_demand_period;
    dtalog.output() << "_____________" << endl;
    dtalog.output() << "Step 1: Reading input data" << endl;
    dtalog.output() << "_____________" << endl;

    dtalog.output() << "Step 1.1: Reading section [demand_period] in setting.csv..." << endl;

    parser_demand_period.IsFirstLineHeader = false;
    if (parser_demand_period.OpenCSVFile("settings.csv", false))
    {
        while (parser_demand_period.ReadRecord_Section())
        {
            if (parser_demand_period.SectionName == "[demand_period]")
            {
                CDemand_Period demand_period;

                if (!parser_demand_period.GetValueByFieldName("demand_period_id", demand_period.demand_period_id))
                    break;

                if (!parser_demand_period.GetValueByFieldName("demand_period", demand_period.demand_period))
                {
                    dtalog.output() << "Error: Field demand_period in file demand_period cannot be read." << endl;
                    g_ProgramStop();
                }

                vector<float> global_minute_vector;

                if (!parser_demand_period.GetValueByFieldName("time_period", demand_period.time_period))
                {
                    dtalog.output() << "Error: Field time_period in file demand_period cannot be read." << endl;
                    g_ProgramStop();
                }

                //input_string includes the start and end time of a time period with hhmm format
                global_minute_vector = g_time_parser(demand_period.time_period); //global_minute_vector incldue the starting and ending time

                if (global_minute_vector.size() == 2)
                {
                    demand_period.starting_time_slot_no = global_minute_vector[0] / MIN_PER_TIMESLOT;
                    demand_period.ending_time_slot_no = global_minute_vector[1] / MIN_PER_TIMESLOT;

                    if (global_minute_vector[0] < assignment.g_LoadingStartTimeInMin)
                        assignment.g_LoadingStartTimeInMin = global_minute_vector[0];

                    if (global_minute_vector[1] > assignment.g_LoadingEndTimeInMin)
                        assignment.g_LoadingEndTimeInMin = global_minute_vector[1];

                    if (assignment.g_LoadingEndTimeInMin < assignment.g_LoadingStartTimeInMin)
                    {
                        assignment.g_LoadingEndTimeInMin = assignment.g_LoadingStartTimeInMin + 1; // in case user errror
                    }
                    //g_fout << global_minute_vector[0] << endl;
                    //g_fout << global_minute_vector[1] << endl;


                    char time_interval_field_name[20];

                    for (int s = max(0, demand_period.starting_time_slot_no - 1); s <= demand_period.ending_time_slot_no; s++)
                    {
                        demand_period.cumulative_departure_time_ratio[s] = 0;
                    }


                    for (int s = demand_period.starting_time_slot_no; s <= demand_period.ending_time_slot_no; s++)
                    {
                        sprintf(time_interval_field_name, "time_interval%d", s);
                        demand_period.departure_time_ratio[s] = 1.0 / (demand_period.ending_time_slot_no - demand_period.starting_time_slot_no + 1);
                        parser_demand_period.GetValueByFieldName(time_interval_field_name, demand_period.departure_time_ratio[s], false);
                    }

                    demand_period.compute_cumulative_profile(demand_period.starting_time_slot_no, demand_period.ending_time_slot_no);

                }

                assignment.demand_period_to_seqno_mapping[demand_period.demand_period] = assignment.g_DemandPeriodVector.size();
                assignment.g_DemandPeriodVector.push_back(demand_period);
            }
        }

        parser_demand_period.CloseCSVFile();

        if (assignment.g_DemandPeriodVector.size() == 0)
        {
            dtalog.output() << "Error:  Section demand_period has no information." << endl;
            g_ProgramStop();
        }
    }
    else
    {
        dtalog.output() << "Error: File settings.csv cannot be opened.\n It might be currently used and locked by EXCEL." << endl;
        g_ProgramStop();
    }

    dtalog.output() << "number of demand periods = " << assignment.g_DemandPeriodVector.size() << endl;

    assignment.g_number_of_demand_periods = assignment.g_DemandPeriodVector.size();
    //step 1:read demand type file

    dtalog.output() << "Step 1.2: Reading section [link_type] in setting.csv..." << endl;

    CCSVParser parser_link_type;
    parser_link_type.IsFirstLineHeader = false;
    if (parser_link_type.OpenCSVFile("settings.csv", false))
    {
        // create a special link type as virtual connector
        CLinkType element_vc;
        // -1 is for virutal connector
        element_vc.link_type = -1;
        element_vc.type_code = "c";
        element_vc.traffic_flow_code = 0;
        assignment.g_LinkTypeMap[element_vc.link_type] = element_vc;
        //end of create special link type for virtual connectors

        int line_no = 0;

        while (parser_link_type.ReadRecord_Section())
        {
            if (parser_link_type.SectionName == "[link_type]")
            {
                CLinkType element;

                if (!parser_link_type.GetValueByFieldName("link_type", element.link_type))
                {
                    if (line_no == 0)
                    {
                        dtalog.output() << "Error: Field link_type cannot be found in file link_type.csv." << endl;
                        g_ProgramStop();
                    }
                    else
                    {
                        // read empty line
                        break;
                    }
                }

                if (assignment.g_LinkTypeMap.find(element.link_type) != assignment.g_LinkTypeMap.end())
                {
                    dtalog.output() << "Error: Field link_type " << element.link_type << " has been defined more than once in file link_type.csv." << endl;
                    g_ProgramStop();
                }

                string traffic_flow_code_str;
                parser_link_type.GetValueByFieldName("type_code", element.type_code, true);
                parser_link_type.GetValueByFieldName("traffic_flow_code", traffic_flow_code_str);

                // by default bpr
                element.traffic_flow_code = 0;

                if (traffic_flow_code_str == "point_queue")
                    element.traffic_flow_code = 1;

                if (traffic_flow_code_str == "spatial_queue")
                    element.traffic_flow_code = 2;

                if (traffic_flow_code_str == "kw")
                    element.traffic_flow_code = 3;

                dtalog.output() << "important: traffic_flow_code on link type " << element.link_type << " is " << element.traffic_flow_code << endl;


                assignment.g_LinkTypeMap[element.link_type] = element;
                line_no++;
            }
        }

        parser_link_type.CloseCSVFile();
    }

    dtalog.output() << "number of link types = " << assignment.g_LinkTypeMap.size() << endl;

    CCSVParser parser_agent_type;
    dtalog.output() << "Step 1.3: Reading section [agent_type] in setting.csv..." << endl;

    parser_agent_type.IsFirstLineHeader = false;
    if (parser_agent_type.OpenCSVFile("settings.csv", false))
    {
        assignment.g_AgentTypeVector.clear();
        while (parser_agent_type.ReadRecord_Section())
        {
            if (parser_agent_type.SectionName == "[agent_type]")
            {
                CAgent_type agent_type;
                agent_type.agent_type_no = assignment.g_AgentTypeVector.size();

                if (!parser_agent_type.GetValueByFieldName("agent_type", agent_type.agent_type))
                    break;

                parser_agent_type.GetValueByFieldName("VOT", agent_type.value_of_time, false, false);

                // scan through the map with different node sum for different paths
                parser_agent_type.GetValueByFieldName("PCE", agent_type.PCE, false, false);
                parser_agent_type.GetValueByFieldName("headway", agent_type.time_headway_in_sec, false, false);
                parser_agent_type.GetValueByFieldName("display_code", agent_type.display_code);
                parser_agent_type.GetValueByFieldName("real_time_information_type", agent_type.real_time_information);

                if (agent_type.agent_type == "vms")  // set the real time information type = 1 for vms class by default
                    agent_type.real_time_information = 1;

                parser_agent_type.GetValueByFieldName("access_node_type", agent_type.access_node_type, false);

                if (agent_type.access_node_type.size() > 0)
                {
                    parser_agent_type.GetValueByFieldName("access_speed", agent_type.access_speed);
                    parser_agent_type.GetValueByFieldName("access_distance_lb", agent_type.access_distance_lb);
                    parser_agent_type.GetValueByFieldName("access_distance_ub", agent_type.access_distance_ub);
                    parser_agent_type.GetValueByFieldName("acecss_link_k", agent_type.acecss_link_k);
                }

                assignment.agent_type_2_seqno_mapping[agent_type.agent_type] = assignment.g_AgentTypeVector.size();

                assignment.g_AgentTypeVector.push_back(agent_type);
                assignment.g_number_of_agent_types = assignment.g_AgentTypeVector.size();
            }
        }
        parser_agent_type.CloseCSVFile();

        if (assignment.g_AgentTypeVector.size() == 0)
            dtalog.output() << "Error: Section agent_type does not contain information." << endl;
    }

    if (assignment.g_AgentTypeVector.size() >= _MAX_AGNETTYPES)
    {
        dtalog.output() << "Error: agent_type = " << assignment.g_AgentTypeVector.size() << " in section agent_type is too large. " << "_MAX_AGNETTYPES = " << _MAX_AGNETTYPES << "Please contact program developers!";
        g_ProgramStop();
    }

    dtalog.output() << "number of agent typess = " << assignment.g_AgentTypeVector.size() << endl;

    assignment.g_number_of_nodes = 0;
    assignment.g_number_of_links = 0;  // initialize  the counter to 0

    int internal_node_seq_no = 0;
    // step 3: read node file

    std::map<int, int> zone_id_to_centriod_node_id_mapping;  // this is an one-to-one mapping
    std::map<int, int> zone_id_mapping;  // this is used to mark if this zone_id has been identified or not
    std::map<int, double> zone_id_x;
    std::map<int, double> zone_id_y;
    std::map<int, int> info_zone_id_mapping;  // this is used to mark if this zone_id has been provided with information or not



    std::map<int, int> zone_id_production;
    std::map<int, int> zone_id_attraction;

    CCSVParser parser;

    int multmodal_activity_node_count = 0;
    dtalog.output() << "Step 1.4: Reading node data in node.csv..." << endl;

    if (parser.OpenCSVFile("node.csv", true))
    {
        while (parser.ReadRecord())  // if this line contains [] mark, then we will also read field headers.
        {
            int node_id;
            if (!parser.GetValueByFieldName("node_id", node_id))
                continue;

            if (assignment.g_node_id_to_seq_no_map.find(node_id) != assignment.g_node_id_to_seq_no_map.end())
            {
                //has been defined
                continue;
            }
            assignment.g_node_id_to_seq_no_map[node_id] = internal_node_seq_no;

            // create a node object
            CNode node;
            node.node_id = node_id;
            node.node_seq_no = internal_node_seq_no;

            int zone_id = -1;

            parser.GetValueByFieldName("node_type", node.node_type, false);

            parser.GetValueByFieldName("zone_id", zone_id);
            if (zone_id >= 1)
            {
                node.is_activity_node = 1;  // from zone
                string str_agent_type;
                parser.GetValueByFieldName("agent_type", str_agent_type, false);

                if (str_agent_type.size() > 0 && assignment.agent_type_2_seqno_mapping.find(str_agent_type) != assignment.agent_type_2_seqno_mapping.end())
                {
                    node.agent_type_str = str_agent_type;
                    node.agent_type_no = assignment.agent_type_2_seqno_mapping[str_agent_type];
                    multmodal_activity_node_count++;
                }
            }

            parser.GetValueByFieldName("x_coord", node.x, true, false);
            parser.GetValueByFieldName("y_coord", node.y, true, false);

            int subarea_id = -1;
            parser.GetValueByFieldName("subarea_id", subarea_id, false);
            node.subarea_id = subarea_id;
            // this is an activity node // we do not allow zone id of zero
            if (zone_id >= 1)
            {
                // for physcial nodes because only centriod can have valid zone_id.
                node.zone_org_id = zone_id;
                if (zone_id_mapping.find(zone_id) == zone_id_mapping.end())
                {
                    //create zone
                    zone_id_mapping[zone_id] = node_id;

                    assignment.zone_id_X_mapping[zone_id] = node.x;
                    assignment.zone_id_Y_mapping[zone_id] = node.y;
                }

                // for od calibration, I think we don't need to implement for now
                if (assignment.assignment_mode == 5)
                {
                    float production = 0;
                    float attraction = 0;
                    parser.GetValueByFieldName("production", production);
                    parser.GetValueByFieldName("attraction", attraction);

                    zone_id_production[zone_id] = production;
                    zone_id_attraction[zone_id] = attraction;
                }
            }
            int info_type = 0;
            parser.GetValueByFieldName("info_type", info_type, false);

            if (info_type >= 1 && zone_id >= 1)
            {

                info_zone_id_mapping[zone_id] = node_id;
                node.is_information_zone = 1;
                assignment.node_seq_no_2_info_zone_id_mapping[internal_node_seq_no] = zone_id;
            }

            /*node.x = x;
            node.y = y;*/
            internal_node_seq_no++;

            // push it to the global node vector
            g_node_vector.push_back(node);
            assignment.g_number_of_nodes++;

            if (assignment.g_number_of_nodes % 5000 == 0)
                dtalog.output() << "reading " << assignment.g_number_of_nodes << " nodes.. " << endl;
        }

        dtalog.output() << "number of nodes = " << assignment.g_number_of_nodes << endl;
        dtalog.output() << "number of multimodal activity nodes = " << multmodal_activity_node_count << endl;
        dtalog.output() << "number of zones = " << zone_id_mapping.size() << endl;
        dtalog.output() << "number of info zones = " << info_zone_id_mapping.size() << endl;

        // fprintf(g_pFileOutputLog, "number of nodes =,%d\n", assignment.g_number_of_nodes);
        parser.CloseCSVFile();
    }


    /// <summary>  mappping node to zone
    // hanlding multimodal access link: stage 1
    for (int at = 0; at < assignment.g_AgentTypeVector.size(); ++at) // first loop for each agent type
    {
        if (assignment.g_AgentTypeVector[at].access_node_type.size() > 0)  // for each multmodal agent type
        {
            // find the closest zone id

            for (int a_k = 0; a_k < g_node_vector.size(); a_k++)
            {
                if (g_node_vector[a_k].is_activity_node == 1 && g_node_vector[a_k].agent_type_no == at) //second loop for mode_specific activity node
                {

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
                        if (g_node_vector[i].node_type.size() > 0)  // stop or station  //third loop for each stop or station node
                        {

                            if (assignment.g_AgentTypeVector[at].access_node_type.find(g_node_vector[i].node_type) != string::npos)  // check allowed access code
                            {

                                double zone_x = g_node_vector[a_k].x;
                                double zone_y = g_node_vector[a_k].y;

                                double distance = g_calculate_p2p_distance_in_mile_from_latitude_longitude(zone_x, zone_y, g_node_vector[i].x, g_node_vector[i].y);
                                // calculate the distance 

                                if (distance < min_distance)
                                {
                                    min_distance = distance;
                                    min_distance_node_seq_no = i;
                                }

                                if (distance >= assignment.g_AgentTypeVector[at].access_distance_lb && distance <= assignment.g_AgentTypeVector[at].access_distance_ub)  // check the range 
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
                        float distance_k = 99999;

                        if(access_node_distance_vector.size() > assignment.g_AgentTypeVector[at].acecss_link_k)
                        {

                        std::vector<float> access_node_distance_vector_temp;
                        access_node_distance_vector_temp = access_node_distance_vector;
                        std::sort(access_node_distance_vector_temp.begin(), access_node_distance_vector_temp.end());

                        distance_k = access_node_distance_vector_temp[max(0,assignment.g_AgentTypeVector[at].acecss_link_k - 1)];

                        }

                        for (int an = 0; an < access_node_seq_vector.size(); an++)
                        {
                            if(access_node_distance_vector[an] < distance_k)  // within the shortest k ranage 
                            { 
                            g_add_new_access_link(a_k, access_node_seq_vector[an], access_node_distance_vector[an], at);
                            g_add_new_access_link(access_node_seq_vector[an], a_k, access_node_distance_vector[an], at);
                            }
                        }

                    }
                    else if (min_distance_node_seq_no >= 0 && min_distance < assignment.g_AgentTypeVector[at].access_distance_ub)  // no node in the  preferred range, just use any feasible node with minimum distance by default
                    {
                        g_add_new_access_link(a_k, min_distance_node_seq_no, min_distance, at);
                        g_add_new_access_link(min_distance_node_seq_no, a_k, min_distance, at);

                    }
                    else {

//                        dtalog.output() << " zone" << g_node_vector[a_k].zone_org_id << " with agent type = " << assignment.g_AgentTypeVector[at].agent_type.c_str() << " has no access to stop or station" << endl;
                        // zone without multimodal access
                    }

                }  // for each zone

            }
        } 
    }// for each agent type 

        //g_InfoZoneMapping(assignment);
        g_OutputModelFiles(1);  // node

        // initialize zone vector
        dtalog.output() << "Step 1.5: Initializing O-D zone vector..." << endl;

        std::map<int, int>::iterator it;
        // creating zone centriod
        for (it = zone_id_mapping.begin(); it != zone_id_mapping.end(); ++it)
        {
            COZone ozone;

            // for each zone, we have to also create centriod
            ozone.zone_id = it->first;  // zone_id
            ozone.zone_seq_no = g_zone_vector.size();
            ozone.obs_production = zone_id_production[it->first];
            ozone.obs_attraction = zone_id_attraction[it->first];

            assignment.g_zoneid_to_zone_seq_no_mapping[ozone.zone_id] = ozone.zone_seq_no;  // create the zone id to zone seq no mapping

            // create a centriod
            CNode node;
            // very large number as a special id
            node.node_id = -1 * ozone.zone_id;
            node.node_seq_no = g_node_vector.size();
            assignment.g_node_id_to_seq_no_map[node.node_id] = node.node_seq_no;
            node.zone_id = ozone.zone_id;

            if (info_zone_id_mapping.find(ozone.zone_id) != info_zone_id_mapping.end())
            {
                ozone.b_real_time_information = true;
                assignment.zone_seq_no_2_info_mapping[ozone.zone_seq_no] = 1;
                dtalog.output() << "info zone =" << ozone.zone_id << endl;

            }
            // push it to the global node vector
            g_node_vector.push_back(node);
            assignment.g_number_of_nodes++;

            ozone.node_seq_no = node.node_seq_no;
            // this should be the only one place that defines this mapping
            zone_id_to_centriod_node_id_mapping[ozone.zone_id] = node.node_id;
            // add element into vector
            g_zone_vector.push_back(ozone);
        }

        // gravity model.
        if (assignment.assignment_mode == 5)
        {
            dtalog.output() << "writing demand.csv.." << endl;

            FILE* g_pFileODMatrix = nullptr;
            fopen_ss(&g_pFileODMatrix, "demand.csv", "w");

            if (!g_pFileODMatrix)
            {
                dtalog.output() << "File demand.csv cannot be opened." << endl;
                g_ProgramStop();
            }
            else
            {
                fprintf(g_pFileODMatrix, "o_zone_id,d_zone_id,volume\n");

                float total_attraction = 0;

                for (int d = 0; d < g_zone_vector.size(); ++d)
                {
                    if (g_zone_vector[d].obs_attraction > 0)
                        total_attraction += g_zone_vector[d].obs_attraction;
                }

                // reset the estimated production and attraction
                for (int orig = 0; orig < g_zone_vector.size(); ++orig)  // o
                {
                    if (g_zone_vector[orig].obs_production >= 0)
                    {
                        for (int dest = 0; dest < g_zone_vector.size(); ++dest)  // d
                        {
                            if (g_zone_vector[dest].obs_attraction > 0)
                            {
                                float value = g_zone_vector[orig].obs_production * g_zone_vector[dest].obs_attraction / max(0.0001f, total_attraction);
                                fprintf(g_pFileODMatrix, "%d,%d,%.4f,\n", g_zone_vector[orig].zone_id, g_zone_vector[dest].zone_id, value);
                                dtalog.output() << "orig= " << g_zone_vector[orig].zone_id << " dest= " << g_zone_vector[dest].zone_id << ":" << value << endl;
                            }
                        }
                    }
                }

                fclose(g_pFileODMatrix);
            }
        }

        dtalog.output() << "number of zones = " << g_zone_vector.size() << endl;
        // step 4: read link file

        CCSVParser parser_link;

        dtalog.output() << "Step 1.6: Reading link data in link.csv... " << endl;
        if (parser_link.OpenCSVFile("link.csv", true))
        {
            while (parser_link.ReadRecord())  // if this line contains [] mark, then we will also read field headers.
            {
                string link_type_name_str;
                parser_link.GetValueByFieldName("link_type_name", link_type_name_str, false);

                if (link_type_name_str == "z2sta"/* || link_type_name_str == "sta2sta_2r"*/)  // filtering
                    continue;

                int from_node_id;
                if (!parser_link.GetValueByFieldName("from_node_id", from_node_id))
                    continue;

                int to_node_id;
                if (!parser_link.GetValueByFieldName("to_node_id", to_node_id))
                    continue;

                string linkID;
                parser_link.GetValueByFieldName("link_id", linkID);
                // add the to node id into the outbound (adjacent) node list

                if (assignment.g_node_id_to_seq_no_map.find(from_node_id) == assignment.g_node_id_to_seq_no_map.end())
                {
                    dtalog.output() << "Error: from_node_id " << from_node_id << " in file link.csv is not defined in node.csv." << endl;
                    continue; //has not been defined
                }

                if (assignment.g_node_id_to_seq_no_map.find(to_node_id) == assignment.g_node_id_to_seq_no_map.end())
                {
                    dtalog.output() << "Error: to_node_id " << to_node_id << " in file link.csv is not defined in node.csv." << endl;
                    continue; //has not been defined
                }

                //if (assignment.g_link_id_map.find(linkID) != assignment.g_link_id_map.end())
                //    dtalog.output() << "Error: link_id " << linkID.c_str() << " has been defined more than once. Please check link.csv." << endl;

                int internal_from_node_seq_no = assignment.g_node_id_to_seq_no_map[from_node_id];  // map external node number to internal node seq no.
                int internal_to_node_seq_no = assignment.g_node_id_to_seq_no_map[to_node_id];


                // create a link object
                CLink link;

                link.from_node_seq_no = internal_from_node_seq_no;
                link.to_node_seq_no = internal_to_node_seq_no;
                link.link_seq_no = assignment.g_number_of_links;
                link.to_node_seq_no = internal_to_node_seq_no;
                link.link_id = linkID;

                assignment.g_link_id_map[link.link_id] = 1;

                string movement_str;
                parser_link.GetValueByFieldName("mvmt_txt_id", movement_str, false);
                int cell_type = -1;
                if (parser_link.GetValueByFieldName("cell_type", cell_type, false) == true)
                    link.cell_type = cell_type;




                parser_link.GetValueByFieldName("geometry", link.geometry, false);
                parser_link.GetValueByFieldName("path_code", link.path_code_str, false);
                parser_link.GetValueByFieldName("tmc_corridor_name", link.tmc_corridor_name, false);
                parser_link.GetValueByFieldName("link_type_name", link.link_type_name, false);

                parser_link.GetValueByFieldName("link_type_code", link.link_type_code, false);

                // and valid
                if (movement_str.size() > 0)
                {
                    int main_node_id = -1;


                    link.mvmt_txt_id = movement_str;
                    link.main_node_id = main_node_id;
                }

                // Peiheng, 05/13/21, if setting.csv does not have corresponding link type or the whole section is missing, set it as 2 (i.e., Major arterial)
                int link_type = 2;
                parser_link.GetValueByFieldName("link_type", link_type, false);

                if (assignment.g_LinkTypeMap.find(link_type) == assignment.g_LinkTypeMap.end())
                {
                    dtalog.output() << "link type " << link_type << " in link.csv is not defined for link " << from_node_id << "->" << to_node_id << " in link_type.csv" << endl;
                    // link.link_type has been taken care by its default constructor
                    //g_ProgramStop();
                }
                else
                {
                    // link type should be defined in settings.csv
                    link.link_type = link_type;
                }

                if (assignment.g_LinkTypeMap[link.link_type].type_code == "c")  // suggestion: we can move "c" as "connector" in allowed_uses
                {
                    if (g_node_vector[internal_from_node_seq_no].zone_org_id >= 0)
                    {
                        int zone_org_id = g_node_vector[internal_from_node_seq_no].zone_org_id;
                        if (assignment.g_zoneid_to_zone_seq_no_mapping.find(zone_org_id) != assignment.g_zoneid_to_zone_seq_no_mapping.end())
                            link.zone_seq_no_for_outgoing_connector = assignment.g_zoneid_to_zone_seq_no_mapping[zone_org_id];
                    }
                }

                double length = 1.0; // km or mile
                double free_speed = 1.0;
                double k_jam = 200;
                double bwtt_speed = 12;  //miles per hour

                double lane_capacity = 1800;
                parser_link.GetValueByFieldName("length", length);
                if (length < 0.007)
                {
                    length = 0.007;  // minimum length
                }
                parser_link.GetValueByFieldName("free_speed", free_speed);

                if (free_speed <= 0.1)
                    free_speed = 60;

                free_speed = max(0.1, free_speed);

                link.free_speed = free_speed;



                int number_of_lanes = 1;
                parser_link.GetValueByFieldName("lanes", number_of_lanes);
                parser_link.GetValueByFieldName("capacity", lane_capacity);

                link.free_flow_travel_time_in_min = length / free_speed * 60;
                link.traffic_flow_code = assignment.g_LinkTypeMap[link.link_type].traffic_flow_code;

                //spatial queue and kinematic wave
                link.spatial_capacity_in_vehicles = max(1.0, length * number_of_lanes * k_jam);

                // kinematic wave
                if (link.traffic_flow_code == 3)
                    link.BWTT_in_simulation_interval = length / bwtt_speed * 3600 / number_of_seconds_per_interval;

                // Peiheng, 02/03/21, useless block
                if (linkID == "10")
                    int i_debug = 1;

                char VDF_field_name[50];

                for (int at = 0; at < assignment.g_AgentTypeVector.size(); at++)
                {
                    double pce_at = 1; // default
                    sprintf(VDF_field_name, "VDF_pce%s", assignment.g_AgentTypeVector[at].agent_type.c_str());

                    parser_link.GetValueByFieldName(VDF_field_name, pce_at, false, true);

                    if (pce_at > 1.001)  // log
                    {
                        //dtalog.output() << "link " << from_node_id << "->" << to_node_id << " has a pce of " << pce_at << " for agent type "
                        //    << assignment.g_AgentTypeVector[at].agent_type.c_str() << endl;
                    }


                    for (int tau = 0; tau < assignment.g_number_of_demand_periods; ++tau)
                    {
                        link.VDF_period[tau].pce[at] = pce_at;
                    }

                }

                // for traffic simulation
                for (int s = 0; s < _MAX_TIMESLOT_PerPeriod; s++)
                {
                    link.dynamic_link_capacity[s] = lane_capacity * number_of_lanes;;
                }

                for (int tau = 0; tau < assignment.g_number_of_demand_periods; ++tau)
                {
                    //setup default values
                    link.VDF_period[tau].capacity = lane_capacity * number_of_lanes;
                    link.VDF_period[tau].FFTT = length / free_speed * 60.0;  // 60.0 for 60 min per hour
                    link.VDF_period[tau].alpha = 0.15;
                    link.VDF_period[tau].beta = 4;
                    link.VDF_period[tau].preload = 0;

                    for (int at = 0; at < assignment.g_AgentTypeVector.size(); at++)
                    {
                        link.VDF_period[tau].toll[at] = 0;
                        link.VDF_period[tau].LR_price[at] = 0;
                        link.VDF_period[tau].LR_RT_price[at] = 0;
                    }


                    link.VDF_period[tau].starting_time_slot_no = assignment.g_DemandPeriodVector[tau].starting_time_slot_no;
                    link.VDF_period[tau].ending_time_slot_no = assignment.g_DemandPeriodVector[tau].ending_time_slot_no;


                    int demand_period_id = assignment.g_DemandPeriodVector[tau].demand_period_id;
                    sprintf(VDF_field_name, "VDF_fftt%d", demand_period_id);
                    parser_link.GetValueByFieldName(VDF_field_name, link.VDF_period[tau].FFTT, false, false);  // FFTT should be per min

                    if (link.VDF_period[tau].FFTT > 100)
                    {
                        dtalog.output() << "link " << from_node_id << "->" << to_node_id << " has a FFTT of " << link.VDF_period[tau].FFTT << " min at demand period " << demand_period_id
                            << " " << assignment.g_DemandPeriodVector[tau].demand_period.c_str() << endl;
                    }

                    sprintf(VDF_field_name, "VDF_cap%d", demand_period_id);
                    parser_link.GetValueByFieldName(VDF_field_name, link.VDF_period[tau].capacity, false, false);  // capacity should be per period per link (include all lanes)

                    sprintf(VDF_field_name, "VDF_alpha%d", demand_period_id);
                    parser_link.GetValueByFieldName(VDF_field_name, link.VDF_period[tau].alpha, false, false);

                    sprintf(VDF_field_name, "VDF_beta%d", demand_period_id);
                    parser_link.GetValueByFieldName(VDF_field_name, link.VDF_period[tau].beta, false, false);

                    sprintf(VDF_field_name, "VDF_allowed_uses%d", demand_period_id);
                    parser_link.GetValueByFieldName(VDF_field_name, link.VDF_period[tau].allowed_uses, false);

                    sprintf(VDF_field_name, "VDF_preload%d", demand_period_id);
                    parser_link.GetValueByFieldName(VDF_field_name, link.VDF_period[tau].preload, false, false);

                    for (int at = 0; at < assignment.g_AgentTypeVector.size(); at++)
                    {
                        sprintf(VDF_field_name, "VDF_toll%s%d", assignment.g_AgentTypeVector[at].agent_type.c_str(), demand_period_id);
                        parser_link.GetValueByFieldName(VDF_field_name, link.VDF_period[tau].toll[at], false, false);

                        if (link.VDF_period[tau].toll[at] > 0.001)
                        {
                            dtalog.output() << "link " << from_node_id << "->" << to_node_id << " has a toll of " << link.VDF_period[tau].toll[at] << " for agent type "
                                << assignment.g_AgentTypeVector[at].agent_type.c_str() << " at demand period " << demand_period_id << endl;
                        }

                    }

                    sprintf(VDF_field_name, "VDF_penalty%d", demand_period_id);
                    parser_link.GetValueByFieldName(VDF_field_name, link.VDF_period[tau].penalty, false, false);

                    if (link.cell_type >= 1) // micro lane-changing arc
                    {
                        // additinonal min: 24 seconds 0.4 min
                        link.VDF_period[tau].penalty += 0.4;
                    }

                    sprintf(VDF_field_name, "VDF_PHF%d", demand_period_id);
                    parser_link.GetValueByFieldName(VDF_field_name, link.VDF_period[tau].PHF, false, false);

                    sprintf(VDF_field_name, "VDF_mu%d", demand_period_id);
                    parser_link.GetValueByFieldName(VDF_field_name, link.VDF_period[tau].mu, false, false);  // mu should be per hour per link, so that we can calculate congestion duration and D/mu in BPR-X

                    parser_link.GetValueByFieldName("cycle_length", link.VDF_period[tau].cycle_length, false, false);

                    if (link.VDF_period[tau].cycle_length >= 1)
                    {
                        link.timing_arc_flag = true;

                        parser_link.GetValueByFieldName("start_green_time", link.VDF_period[tau].start_green_time);
                        parser_link.GetValueByFieldName("end_green_time", link.VDF_period[tau].end_green_time);
                        parser_link.GetValueByFieldName("red_time", link.VDF_period[tau].red_time, false);
                        parser_link.GetValueByFieldName("green_time", link.VDF_period[tau].effective_green_time, false);
                    }

                }

                // for each period

                float default_cap = 1000;
                float default_BaseTT = 1;

                // setup default value
                for (int tau = 0; tau < assignment.g_number_of_demand_periods; ++tau)
                {
                    link.TDBaseTT[tau] = default_BaseTT;
                    link.TDBaseCap[tau] = default_cap;
                }

                //link.m_OutflowNumLanes = number_of_lanes;//visum lane_cap is actually link_cap

                link.number_of_lanes = number_of_lanes;
                link.lane_capacity = lane_capacity;
                link.link_spatial_capacity = k_jam * number_of_lanes * length;

                link.length = max(0.00001, length);
                for (int tau = 0; tau < assignment.g_number_of_demand_periods; ++tau)
                    link.travel_time_per_period[tau] = length / free_speed * 60;

                // min // calculate link cost based length and speed limit // later we should also read link_capacity, calculate delay

                //int sequential_copying = 0;
                //
                //parser_link.GetValueByFieldName("sequential_copying", sequential_copying);

                g_node_vector[internal_from_node_seq_no].m_outgoing_link_seq_no_vector.push_back(link.link_seq_no);  // add this link to the corresponding node as part of outgoing node/link
                g_node_vector[internal_to_node_seq_no].m_incoming_link_seq_no_vector.push_back(link.link_seq_no);  // add this link to the corresponding node as part of outgoing node/link

                g_node_vector[internal_from_node_seq_no].m_to_node_seq_no_vector.push_back(link.to_node_seq_no);  // add this link to the corresponding node as part of outgoing node/link
                g_node_vector[internal_from_node_seq_no].m_to_node_2_link_seq_no_map[link.to_node_seq_no] = link.link_seq_no;  // add this link to the corresponding node as part of outgoing node/link

                g_link_vector.push_back(link);

                string mvmt_key;
                parser_link.GetValueByFieldName("mvmt_key", mvmt_key, false);
                if (mvmt_key.size() > 4) // main_node_id _ movement code
                {
                    assignment.g_mvmt_key_to_link_no_map[mvmt_key] = assignment.g_number_of_links;
                }

                assignment.g_number_of_links++;

                if (assignment.g_number_of_links % 10000 == 0)
                    dtalog.output() << "reading " << assignment.g_number_of_links << " links.. " << endl;
            }

            parser_link.CloseCSVFile();
        }
        // we now know the number of links
        dtalog.output() << "number of links = " << assignment.g_number_of_links << endl;

        // after we read the physical links
        // we create virtual connectors
        for (int i = 0; i < g_node_vector.size(); ++i)
        {

            if (g_node_vector[i].zone_org_id >= 0) // for each physical node
            { // we need to make sure we only create two way connectors between nodes and zones

                    // for each node-zone pair: create a pair of connectors with the agent-type related acess_map
                int zone_org_id = g_node_vector[i].zone_org_id;
                int internal_from_node_seq_no, internal_to_node_seq_no, zone_seq_no;

                internal_from_node_seq_no = g_node_vector[i].node_seq_no;
                int node_id = zone_id_to_centriod_node_id_mapping[zone_org_id];
                internal_to_node_seq_no = assignment.g_node_id_to_seq_no_map[node_id];
                zone_seq_no = assignment.g_zoneid_to_zone_seq_no_mapping[zone_org_id];

                // we need to mark all accessble model on this access links, so we can handle that in the future for each agent type's memory block in shortest path
                // incomming virtual connector
                g_add_new_virtual_connector_link(internal_from_node_seq_no, internal_to_node_seq_no, g_node_vector[i].agent_type_str, -1);
                // outgoing virtual connector
                g_add_new_virtual_connector_link(internal_to_node_seq_no, internal_from_node_seq_no, g_node_vector[i].agent_type_str, zone_seq_no);
                // result is that: we have a unique pair of node-zone access link in the overall network, but still carry agent_type_acess_map for agent types with access on this node-zone connector

            }

        }
        dtalog.output() << "number of links =" << assignment.g_number_of_links << endl;
        g_OutputModelFiles(3);

        if (dtalog.debug_level() == 2)
        {
            for (int i = 0; i < g_node_vector.size(); ++i)
            {
                if (g_node_vector[i].zone_org_id > 0) // for each physical node
                {
                    // we need to make sure we only create two way connectors between nodes and zones
                    dtalog.output() << "node id= " << g_node_vector[i].node_id << " with zone id " << g_node_vector[i].zone_org_id << "and "
                        << g_node_vector[i].m_outgoing_link_seq_no_vector.size() << " outgoing links." << endl;
                    for (int j = 0; j < g_node_vector[i].m_outgoing_link_seq_no_vector.size(); ++j)
                    {
                        int link_seq_no = g_node_vector[i].m_outgoing_link_seq_no_vector[j];
                        dtalog.output() << "  outgoing node = " << g_node_vector[g_link_vector[link_seq_no].to_node_seq_no].node_id << endl;
                    }
                }
                else
                {
                    if (dtalog.debug_level() == 3)
                    {
                        dtalog.output() << "node id= " << g_node_vector[i].node_id << " with " << g_node_vector[i].m_outgoing_link_seq_no_vector.size() << " outgoing links." << endl;
                        for (int j = 0; j < g_node_vector[i].m_outgoing_link_seq_no_vector.size(); ++j)
                        {
                            int link_seq_no = g_node_vector[i].m_outgoing_link_seq_no_vector[j];
                            dtalog.output() << "  outgoing node = " << g_node_vector[g_link_vector[link_seq_no].to_node_seq_no].node_id << endl;
                        }
                    }
                }
            }
        }

}
    //CCSVParser parser_movement;
    //int prohibited_count = 0;

    //if (parser_movement.OpenCSVFile("movement.csv", false))  // not required
    //{
    //    while (parser_movement.ReadRecord())
    //    {
    //        string ib_link_id;
    //        int node_id = 0;
    //        string ob_link_id;
    //        int prohibited_flag = 0;

    //        if (!parser_movement.GetValueByFieldName("node_id", node_id))
    //            break;

    //        if (assignment.g_node_id_to_seq_no_map.find(node_id) == assignment.g_node_id_to_seq_no_map.end())
    //        {
    //            dtalog.output() << "Error: node_id " << node_id << " in file movement.csv is not defined in node.csv." << endl;
    //            //has not been defined
    //            continue;
    //        }

    //        parser_movement.GetValueByFieldName("ib_link_id", ib_link_id);
    //        parser_movement.GetValueByFieldName("ob_link_id", ob_link_id);

    //        if (assignment.g_link_id_map.find(ib_link_id) != assignment.g_link_id_map.end())
    //            dtalog.output() << "Error: ib_link_id " << ib_link_id.c_str() << " has not been defined in movement.csv. Please check link.csv." << endl;

    //        if (assignment.g_link_id_map.find(ob_link_id) != assignment.g_link_id_map.end())
    //            dtalog.output() << "Error: ob_link_id " << ob_link_id.c_str() << " has not been defined in movement.csv. Please check link.csv." << endl;

    //        float penalty = 0;
    //        parser_movement.GetValueByFieldName("penalty", penalty);

    //        if (penalty >= 99)
    //        {
    //            string	movement_string;
    //            movement_string = ib_link_id + "->" + ob_link_id;

    //            int node_no = assignment.g_node_id_to_seq_no_map[node_id];
    //            g_node_vector[node_no].prohibited_movement_size++;
    //            g_node_vector[node_no].m_prohibited_movement_string_map[movement_string] = 1;

    //            prohibited_count++;
    //        }
    //    }

    //    dtalog.output() << "Step XX: Reading movement.csv data with " << prohibited_count << " prohibited records." << endl;
    //    parser_movement.CloseCSVFile();
    //}



//
//
//void g_reload_timing_arc_data(Assignment& assignment)
//{
//    dtalog.output() << "Step 1.7: Reading service arc in timing.csv..." << endl;
//
//    CCSVParser parser_timing_arc;
//    if (parser_timing_arc.OpenCSVFile("timing.csv", false))
//    {
//        while (parser_timing_arc.ReadRecord())  // if this line contains [] mark, then we will also read field headers.
//        {
//            string mvmt_key;
//            if (!parser_timing_arc.GetValueByFieldName("mvmt_key", mvmt_key))
//            {
//                dtalog.output() << "Error: mvmt_key in file timing.csv is not defined." << endl;
//                continue;
//            }
//            // create a link object
//            CSignalTiming timing_arc;
//
//            if (assignment.g_mvmt_key_to_link_no_map.find(mvmt_key) == assignment.g_mvmt_key_to_link_no_map.end())
//            {
//                dtalog.output() << "Error: mvmt_key " << mvmt_key << " in file timing.csv is not defined in link.csv." << endl;
//                //has not been defined
//                continue;
//            }
//            else
//            {
//                timing_arc.link_seq_no = assignment.g_mvmt_key_to_link_no_map[mvmt_key];
//                g_link_vector[timing_arc.link_seq_no].timing_arc_flag = true;
//            }
//
//            string time_period;
//            if (!parser_timing_arc.GetValueByFieldName("time_window", time_period))
//            {
//                dtalog.output() << "Error: Field time_window in file timing.csv cannot be read." << endl;
//                g_ProgramStop();
//                break;
//            }
//
//            vector<float> global_minute_vector;
//
//            //input_string includes the start and end time of a time period with hhmm format
//            global_minute_vector = g_time_parser(time_period); //global_minute_vector incldue the starting and ending time
//            if (global_minute_vector.size() == 2)
//            {
//                if (global_minute_vector[0] < assignment.g_LoadingStartTimeInMin)
//                    global_minute_vector[0] = assignment.g_LoadingStartTimeInMin;
//
//                if (global_minute_vector[0] > assignment.g_LoadingEndTimeInMin)
//                    global_minute_vector[0] = assignment.g_LoadingEndTimeInMin;
//
//                if (global_minute_vector[1] < assignment.g_LoadingStartTimeInMin)
//                    global_minute_vector[1] = assignment.g_LoadingStartTimeInMin;
//
//                if (global_minute_vector[1] > assignment.g_LoadingEndTimeInMin)
//                    global_minute_vector[1] = assignment.g_LoadingEndTimeInMin;
//
//                if (global_minute_vector[1] < global_minute_vector[0])
//                    global_minute_vector[1] = global_minute_vector[0];
//
//            }
//            else
//                continue;
//
//            float time_interval = 0;
//
//
//            // capacity in the space time arcs
//            float capacity = 1;
//            parser_timing_arc.GetValueByFieldName("capacity", capacity);
//            timing_arc.VDF_capacity = max(0.0f, capacity);
//
//            // capacity in the space time arcs
//            parser_timing_arc.GetValueByFieldName("cycle_length", timing_arc.cycle_length);
//
//            // capacity in the space time arcs
//            parser_timing_arc.GetValueByFieldName("red_time", timing_arc.red_time);
//            parser_timing_arc.GetValueByFieldName("start_green_time", timing_arc.start_green_time);
//            parser_timing_arc.GetValueByFieldName("end_green_time", timing_arc.end_green_time);
//
//            for (int tau = 0; tau < assignment.g_number_of_demand_periods; ++tau)
//            {
//                    // to do: we need to consider multiple periods in the future, Xuesong Zhou, August 20, 2020.
//                g_link_vector[timing_arc.link_seq_no].VDF_period[tau].red_time = timing_arc.red_time;
//                g_link_vector[timing_arc.link_seq_no].VDF_period[tau].cycle_length = timing_arc.cycle_length;
//            }
//
//
//            g_signal_timing_arc_vector.push_back(timing_arc);
//            assignment.g_number_of_timing_arcs++;
//
//            if (assignment.g_number_of_timing_arcs % 10000 == 0)
//                dtalog.output() << "reading " << assignment.g_number_of_timing_arcs << " timing_arcs.. " << endl;
//        }
//
//        parser_timing_arc.CloseCSVFile();
//    }
//
//    dtalog.output() << endl;
//    dtalog.output() << "Step 1.8: Reading file section [demand_file_list] in setting.csv..." << endl;
//    // we now know the number of links
//    dtalog.output() << "number of timing records = " << assignment.g_number_of_timing_arcs << endl << endl;
//}
//

