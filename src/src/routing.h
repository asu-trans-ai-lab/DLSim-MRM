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


struct CNodeForwardStar {
    CNodeForwardStar() : OutgoingLinkNoArray{ nullptr }, OutgoingNodeNoArray{ nullptr }, OutgoingLinkSize{ 0 }
    {
    }

    // Peiheng, 03/22/21, release memory
    ~CNodeForwardStar()
    {
        delete[] OutgoingLinkNoArray;
        delete[] OutgoingNodeNoArray;
    }

    int* OutgoingLinkNoArray;
    int* OutgoingNodeNoArray;
    int OutgoingLinkSize;
};

class NetworkForSP  // mainly for shortest path calculation
{
public:
    NetworkForSP() : temp_path_node_vector_size{ _MAX_LINK_SIZE_IN_A_PATH }, m_value_of_time{ 10 }, bBuildNetwork{ false }, m_memory_block_no{ 0 }, m_agent_type_no{ 0 }, m_tau{ 0 }, b_real_time_information{ false }
    {
    }

    ~NetworkForSP()
    {
        if (m_SENodeList)  //1
            delete[] m_SENodeList;

        if (m_node_status_array)  //2
            delete[] m_node_status_array;

        if (m_label_time_array)  //3
            delete[] m_label_time_array;

        if (m_label_distance_array)  //4
            delete[] m_label_distance_array;

        if (m_node_predecessor)  //5
            delete[] m_node_predecessor;

        if (m_link_predecessor)  //6
            delete[] m_link_predecessor;

        if (m_node_label_cost)  //7
            delete[] m_node_label_cost;

        if (m_link_flow_volume_array)  //8
            delete[] m_link_flow_volume_array;

        if (m_link_genalized_cost_array) //9
            delete[] m_link_genalized_cost_array;

        if (m_link_outgoing_connector_zone_seq_no_array) //10
            delete[] m_link_outgoing_connector_zone_seq_no_array;

        // Peiheng, 03/22/21, memory release on OutgoingLinkNoArray and OutgoingNodeNoArray
        // is taken care by the destructor of CNodeForwardStar
        if (NodeForwardStarArray)
            delete[] NodeForwardStarArray;
    }

    int temp_path_node_vector_size;
    float m_value_of_time;
    bool bBuildNetwork;
    int m_memory_block_no;

    //node seq vector for each ODK
    int temp_path_node_vector[_MAX_LINK_SIZE_IN_A_PATH];
    //node seq vector for each ODK
    int temp_path_link_vector[_MAX_LINK_SIZE_IN_A_PATH];

    bool m_bSingleSP_Flag;

    // assigned nodes for computing
    std::vector<int> m_origin_node_vector;
    std::vector<int> m_origin_zone_seq_no_vector;

    bool b_real_time_information;
    int m_tau; // assigned nodes for computing
    int m_agent_type_no; // assigned nodes for computing

    CNodeForwardStar* NodeForwardStarArray;

    int m_threadNo;  // internal thread number

    int m_ListFront; // used in coding SEL
    int m_ListTail;  // used in coding SEL
    int* m_SENodeList; // used in coding SEL

    // label cost for shortest path calcuating
    double* m_node_label_cost;
    // time-based cost
    double* m_label_time_array;
    // distance-based cost
    double* m_label_distance_array;

    // predecessor for nodes
    int* m_node_predecessor;
    // update status
    int* m_node_status_array;
    // predecessor for this node points to the previous link that updates its label cost (as part of optimality condition) (for easy referencing)
    int* m_link_predecessor;

    double* m_link_flow_volume_array;

    double* m_link_genalized_cost_array;
    int* m_link_outgoing_connector_zone_seq_no_array;

    // major function 1:  allocate memory and initialize the data
    void AllocateMemory(int number_of_nodes, int number_of_links)
    {
        NodeForwardStarArray = new CNodeForwardStar[number_of_nodes];

        m_SENodeList = new int[number_of_nodes];  //1

        m_LinkBasedSEList = new int[number_of_links];  //1;  // dimension: number of links

        m_node_status_array = new int[number_of_nodes];  //2
        m_label_time_array = new double[number_of_nodes];  //3
        m_label_distance_array = new double[number_of_nodes];  //4
        m_node_predecessor = new int[number_of_nodes];  //5
        m_link_predecessor = new int[number_of_nodes];  //6
        m_node_label_cost = new double[number_of_nodes];  //7

        m_link_flow_volume_array = new double[number_of_links];  //8

        m_link_genalized_cost_array = new double[number_of_links];  //9
        m_link_outgoing_connector_zone_seq_no_array = new int[number_of_links]; //10
    }

    bool UpdateGeneralizedLinkCost(int agent_type_no , Assignment* p_assignment, int origin_zone)
    {
        bool negative_cost_flag = false;
        for (int i = 0; i < g_link_vector.size(); ++i)
        {
            CLink* pLink = &(g_link_vector[i]);

            double LR_price = pLink->VDF_period[m_tau].LR_price[agent_type_no];

            if (p_assignment->g_AgentTypeVector[agent_type_no].real_time_information >= 1
                && p_assignment->zone_seq_no_2_info_mapping.find(origin_zone) != p_assignment->zone_seq_no_2_info_mapping.end()
                ) // for agent type with information and the origin zone is information zone
            {
                LR_price = pLink->VDF_period[m_tau].LR_RT_price[agent_type_no];

                if (LR_price < -0.01)
                    negative_cost_flag = true;

                if (fabs(LR_price) > 0.001)
                {
                    dtalog.output() << "link " << pLink->link_id.c_str() << " has a lr RT price of " << LR_price << " for agent type "
                        << assignment.g_AgentTypeVector[agent_type_no].agent_type.c_str() << " at demand period =" << m_tau << endl;
                }

            }

             m_link_genalized_cost_array[i] = pLink->travel_time_per_period[m_tau] + pLink->VDF_period[m_tau].penalty + 
                 LR_price +
                 pLink->VDF_period[m_tau].toll[agent_type_no] / m_value_of_time * 60 + pLink->RT_travel_time;  // *60 as 60 min per hour
            //route_choice_cost 's unit is min
        }

        return negative_cost_flag;
    }

    void BuildNetwork(Assignment* p_assignment)
    {
        if (bBuildNetwork)
            return;

        int m_outgoing_link_seq_no_vector[_MAX_LINK_SIZE_FOR_A_NODE];
        int m_to_node_seq_no_vector[_MAX_LINK_SIZE_FOR_A_NODE];

        for (int i = 0; i < g_link_vector.size(); ++i)
        {
            CLink* pLink = &(g_link_vector[i]);
            m_link_outgoing_connector_zone_seq_no_array[i] = pLink->zone_seq_no_for_outgoing_connector;
        }

        for (int i = 0; i < p_assignment->g_number_of_nodes; ++i) //Initialization for all non-origin nodes
        {
            int outgoing_link_size = 0;

            for (int j = 0; j < g_node_vector[i].m_outgoing_link_seq_no_vector.size(); ++j)
            {

                int link_seq_no = g_node_vector[i].m_outgoing_link_seq_no_vector[j];
                // only predefined allowed agent type can be considered
                if (g_link_vector[link_seq_no].AllowAgentType(p_assignment->g_AgentTypeVector[m_agent_type_no].agent_type, m_tau))
                {
                    m_outgoing_link_seq_no_vector[outgoing_link_size] = link_seq_no;
                    m_to_node_seq_no_vector[outgoing_link_size] = g_node_vector[i].m_to_node_seq_no_vector[j];

                    outgoing_link_size++;

                    if (outgoing_link_size >= _MAX_LINK_SIZE_FOR_A_NODE)
                    {
                        dtalog.output() << " Error: outgoing_link_size >= _MAX_LINK_SIZE_FOR_A_NODE" << endl;
                        // output the log

                        g_OutputModelFiles(1);

                        g_ProgramStop();
                    }
                }
            }

            int node_seq_no = g_node_vector[i].node_seq_no;
            NodeForwardStarArray[node_seq_no].OutgoingLinkSize = outgoing_link_size;

            if (outgoing_link_size >= 1)
            {
                NodeForwardStarArray[node_seq_no].OutgoingLinkNoArray = new int[outgoing_link_size];
                NodeForwardStarArray[node_seq_no].OutgoingNodeNoArray = new int[outgoing_link_size];
            }

            for (int j = 0; j < outgoing_link_size; ++j)
            {
                NodeForwardStarArray[node_seq_no].OutgoingLinkNoArray[j] = m_outgoing_link_seq_no_vector[j];
                NodeForwardStarArray[node_seq_no].OutgoingNodeNoArray[j] = m_to_node_seq_no_vector[j];
            }
        }

        // after dynamic arrays are created for forward star
        if (dtalog.debug_level() == 2)
        {
            dtalog.output() << "add outgoing link data into dynamic array" << endl;

            for (int i = 0; i < g_node_vector.size(); ++i)
            {
                if (g_node_vector[i].zone_org_id > 0) // for each physical node
                { // we need to make sure we only create two way connectors between nodes and zones
                    dtalog.output() << "node id= " << g_node_vector[i].node_id << " with zone id " << g_node_vector[i].zone_org_id << "and "
                        << NodeForwardStarArray[i].OutgoingLinkSize << " outgoing links." << endl;

                    for (int j = 0; j < NodeForwardStarArray[i].OutgoingLinkSize; j++)
                    {
                        int link_seq_no = NodeForwardStarArray[i].OutgoingLinkNoArray[j];
                        dtalog.output() << "  outgoing node = " << g_node_vector[g_link_vector[link_seq_no].to_node_seq_no].node_id << endl;
                    }
                }
                else
                {
                    if (dtalog.debug_level() == 3)
                    {
                        dtalog.output() << "node id= " << g_node_vector[i].node_id << " with "
                            << NodeForwardStarArray[i].OutgoingLinkSize << " outgoing links." << endl;

                        for (int j = 0; j < NodeForwardStarArray[i].OutgoingLinkSize; ++j)
                        {
                            int link_seq_no = NodeForwardStarArray[i].OutgoingLinkNoArray[j];
                            dtalog.output() << "  outgoing node = " << g_node_vector[g_link_vector[link_seq_no].to_node_seq_no].node_id << endl;
                        }
                    }
                }
            }
        }

        m_value_of_time = p_assignment->g_AgentTypeVector[m_agent_type_no].value_of_time;
        bBuildNetwork = true;
    }

    // SEList: scan eligible List implementation: the reason for not using STL-like template is to avoid overhead associated pointer allocation/deallocation
    inline void SEList_clear()
    {
        m_ListFront = -1;
        m_ListTail = -1;
    }

    inline void SEList_push_front(int node)
    {
        if (m_ListFront == -1)  // start from empty
        {
            m_SENodeList[node] = -1;
            m_ListFront = node;
            m_ListTail = node;
        }
        else
        {
            m_SENodeList[node] = m_ListFront;
            m_ListFront = node;
        }
    }

    inline void SEList_push_back(int node)
    {
        if (m_ListFront == -1)  // start from empty
        {
            m_ListFront = node;
            m_ListTail = node;
            m_SENodeList[node] = -1;
        }
        else
        {
            m_SENodeList[m_ListTail] = node;
            m_SENodeList[node] = -1;
            m_ListTail = node;
        }
    }

    inline bool SEList_empty()
    {
        return(m_ListFront == -1);
    }

    //	inline int SEList_front()
    //	{
    //		return m_ListFront;
    //	}

    //	inline void SEList_pop_front()
    //	{
    //		int tempFront = m_ListFront;
    //		m_ListFront = m_SENodeList[m_ListFront];
    //		m_SENodeList[tempFront] = -1;
    //	}

    //major function: update the cost for each node at each SP tree, using a stack from the origin structure

    double backtrace_shortest_path_tree(Assignment& assignment, int iteration_number, int o_node_index);

    //major function 2: // time-dependent label correcting algorithm with double queue implementation
    float optimal_label_correcting(int processor_id, Assignment* p_assignment, int iteration_k, int o_node_index, int d_node_no = -1, bool pure_travel_time_cost = false)
    {
        int local_debugging_flag = 0;
        int SE_loop_count = 0;

        if (iteration_k == 0)
            BuildNetwork(p_assignment);  // based on agent type and link type

        int agent_type = m_agent_type_no; // assigned nodes for computing
        int origin_node = m_origin_node_vector[o_node_index]; // assigned nodes for computing
        int origin_zone = m_origin_zone_seq_no_vector[o_node_index]; // assigned nodes for computing

        bool negative_cost_flag = UpdateGeneralizedLinkCost(agent_type, p_assignment, origin_zone);


        if (negative_cost_flag == true)
        {
            dtalog.output() << "Negative Cost: SP iteration k =  " << iteration_k << ": origin node: " << g_node_vector[origin_node].node_id << endl;
            local_debugging_flag = 1;
            negative_cost_label_correcting(processor_id, p_assignment, iteration_k, o_node_index, d_node_no);
            return true;
        }

        if (p_assignment->g_number_of_nodes >= 1000 && origin_zone % 97 == 0)
            dtalog.output() << "label correcting for zone " << origin_zone << " in processor " << processor_id << endl;

        if (dtalog.debug_level() >= 2)
            dtalog.output() << "SP iteration k =  " << iteration_k << ": origin node: " << g_node_vector[origin_node].node_id << endl;

        int number_of_nodes = p_assignment->g_number_of_nodes;
        //Initialization for all non-origin nodes
        for (int i = 0; i < number_of_nodes; ++i)
        {
            // not scanned
            m_node_status_array[i] = 0;
            m_node_label_cost[i] = _MAX_LABEL_COST;
            // pointer to previous NODE INDEX from the current label at current node and time
            m_link_predecessor[i] = -1;
            // pointer to previous NODE INDEX from the current label at current node and time
            m_node_predecessor[i] = -1;
            // comment out to speed up comuting
            ////m_label_time_array[i] = 0;
            ////m_label_distance_array[i] = 0;
        }

        // int internal_debug_flag = 0;
        if (NodeForwardStarArray[origin_node].OutgoingLinkSize == 0)
            return 0;

        //Initialization for origin node at the preferred departure time, at departure time, cost = 0, otherwise, the delay at origin node
        m_label_time_array[origin_node] = 0;
        m_node_label_cost[origin_node] = 0.0;
        //Mark:	m_label_distance_array[origin_node] = 0.0;

        // Peiheng, 02/05/21, duplicate initialization, remove them later
        // pointer to previous NODE INDEX from the current label at current node and time
        m_link_predecessor[origin_node] = -1;
        // pointer to previous NODE INDEX from the current label at current node and time
        m_node_predecessor[origin_node] = -1;

        SEList_clear();
        SEList_push_back(origin_node);

        int from_node, to_node;
        int link_sqe_no;
        double new_time = 0;
        double new_distance = 0;
        double new_to_node_cost = 0;
        int tempFront;
        while (!(m_ListFront == -1))   //SEList_empty()
        {
            // from_node = SEList_front();
            // SEList_pop_front();  // remove current node FromID from the SE list

            from_node = m_ListFront;//pop a node FromID for scanning
            tempFront = m_ListFront;
            m_ListFront = m_SENodeList[m_ListFront];
            m_SENodeList[tempFront] = -1;

            m_node_status_array[from_node] = 2;

            if (dtalog.log_path() >= 2 || local_debugging_flag )
            {
                dtalog.output() << "SP:scan SE node: " << g_node_vector[from_node].node_id << " with "
                    << NodeForwardStarArray[from_node].OutgoingLinkSize << " outgoing link(s). " << endl;
            }
            //scan all outbound nodes of the current node

            int pred_link_seq_no = m_link_predecessor[from_node];

            // for each link (i,j) belong A(i)
            for (int i = 0; i < NodeForwardStarArray[from_node].OutgoingLinkSize; ++i)
            {
                to_node = NodeForwardStarArray[from_node].OutgoingNodeNoArray[i];
                link_sqe_no = NodeForwardStarArray[from_node].OutgoingLinkNoArray[i];

                if (dtalog.log_path() >= 2 || local_debugging_flag)
                    dtalog.output() << "SP:  checking outgoing node " << g_node_vector[to_node].node_id << endl;

                // if(map (pred_link_seq_no, link_sqe_no) is prohibitted )
                //     then continue; //skip this is not an exact solution algorithm for movement

                if (g_node_vector[from_node].prohibited_movement_size >= 1)
                {
                    if (pred_link_seq_no >= 0)
                    {
                        string	movement_string;
                        string ib_link_id = g_link_vector[pred_link_seq_no].link_id;
                        string ob_link_id = g_link_vector[link_sqe_no].link_id;
                        movement_string = ib_link_id + "->" + ob_link_id;

                        if (g_node_vector[from_node].m_prohibited_movement_string_map.find(movement_string) != g_node_vector[from_node].m_prohibited_movement_string_map.end())
                        {
                            dtalog.output() << "prohibited movement " << movement_string << " will not be used " << endl;
                            continue;
                        }
                    }
                }

                //remark: the more complicated implementation can be found in paper Shortest Path Algorithms In Transportation Models: Classical and Innovative Aspects
                //	A note on least time path computation considering delays and prohibitions for intersection movements

                if (m_link_outgoing_connector_zone_seq_no_array[link_sqe_no] >= 0)
                {
                    if (m_link_outgoing_connector_zone_seq_no_array[link_sqe_no] != origin_zone)
                    {
                        // filter out for an outgoing connector with a centriod zone id different from the origin zone seq no
                        continue;
                    }
                }

                //very important: only origin zone can access the outbound connectors,
                //the other zones do not have access to the outbound connectors

                // Mark				new_time = m_label_time_array[from_node] + pLink->travel_time_per_period[tau];
                // Mark				new_distance = m_label_distance_array[from_node] + pLink->length;
                float additional_cost = 0;

                if (g_link_vector[link_sqe_no].RT_travel_time > 1)  // used in real time routing only
                {
                    additional_cost = g_link_vector[link_sqe_no].RT_travel_time;

                    //if (g_link_vector[link_sqe_no].RT_travel_time > 999)
                    //    continue; //skip this link due to closure
                }


                new_to_node_cost = m_node_label_cost[from_node] + m_link_genalized_cost_array[link_sqe_no] + additional_cost;

                if (dtalog.log_path() || local_debugging_flag)
                {
                    dtalog.output() << "SP:  checking from node " << g_node_vector[from_node].node_id
                        << "  to node " << g_node_vector[to_node].node_id << " cost = " << new_to_node_cost << endl;
                }

                if (new_to_node_cost < m_node_label_cost[to_node]) // we only compare cost at the downstream node ToID at the new arrival time t
                {
                    if (dtalog.log_path() || local_debugging_flag)
                    {
                        dtalog.output() << "SP:  updating node: " << g_node_vector[to_node].node_id << " current cost:" << m_node_label_cost[to_node]
                            << " new cost " << new_to_node_cost << endl;
                    }

                    // update cost label and node/time predecessor
                    // m_label_time_array[to_node] = new_time;
                    // m_label_distance_array[to_node] = new_distance;
                    m_node_label_cost[to_node] = new_to_node_cost;
                    // pointer to previous physical NODE INDEX from the current label at current node and time
                    m_node_predecessor[to_node] = from_node;
                    // pointer to previous physical NODE INDEX from the current label at current node and time
                    m_link_predecessor[to_node] = link_sqe_no;

                    if (dtalog.log_path() || local_debugging_flag)
                    {
                        dtalog.output() << "SP: add node " << g_node_vector[to_node].node_id << " new cost:" << new_to_node_cost
                            << " into SE List " << g_node_vector[to_node].node_id << endl;
                    }

                    // deque updating rule for m_node_status_array
                    if (m_node_status_array[to_node] == 0)
                    {
                        ///// SEList_push_back(to_node);
                        ///// begin of inline block
                        if (m_ListFront == -1)  // start from empty
                        {
                            m_ListFront = to_node;
                            m_ListTail = to_node;
                            m_SENodeList[to_node] = -1;
                        }
                        else
                        {
                            m_SENodeList[m_ListTail] = to_node;
                            m_SENodeList[to_node] = -1;
                            m_ListTail = to_node;
                        }
                        ///// end of inline block

                        m_node_status_array[to_node] = 1;
                    }

                    if (m_node_status_array[to_node] == 2)
                    {
                        /////SEList_push_front(to_node);
                        ///// begin of inline block
                        if (m_ListFront == -1)  // start from empty
                        {
                            m_SENodeList[to_node] = -1;
                            m_ListFront = to_node;
                            m_ListTail = to_node;
                        }
                        else
                        {
                            m_SENodeList[to_node] = m_ListFront;
                            m_ListFront = to_node;
                        }
                        ///// end of inline block

                        m_node_status_array[to_node] = 1;
                    }
                }
            }
        }

        if (dtalog.log_path() || local_debugging_flag)
        {
            dtalog.output() << "SPtree at iteration k = " << iteration_k << " origin node: "
                << g_node_vector[origin_node].node_id << endl;

            //Initialization for all non-origin nodes
            for (int i = 0; i < p_assignment->g_number_of_nodes; ++i)
            {
                int node_pred_id = -1;
                int node_pred_no = m_node_predecessor[i];

                if (node_pred_no >= 0)
                    node_pred_id = g_node_vector[node_pred_no].node_id;

                if (m_node_label_cost[i] < 9999)
                {
                    dtalog.output() << "SP node: " << g_node_vector[i].node_id << " label cost " << m_node_label_cost[i] << "time "
                        << m_label_time_array[i] << "node_pred_id " << node_pred_id << endl;
                }
            }
        }

        if (d_node_no >= 1)
            return m_node_label_cost[d_node_no];
        else
            return 0;  // one to all shortest pat
    }

    float negative_cost_label_correcting(int processor_id, Assignment* p_assignment, int iteration_k, int o_node_index, int d_node_no = -1, bool pure_travel_time_cost = false)
    {
        int local_debugging_flag = 0;
        int SE_loop_count = 0;

       int agent_type = m_agent_type_no; // assigned nodes for computing
        int origin_node = m_origin_node_vector[o_node_index]; // assigned nodes for computing
        int origin_zone = m_origin_zone_seq_no_vector[o_node_index]; // assigned nodes for computing

        if (dtalog.debug_level() >= 2)
            dtalog.output() << "SP iteration k =  " << iteration_k << ": origin node: " << g_node_vector[origin_node].node_id << endl;

        int number_of_nodes = p_assignment->g_number_of_nodes;
        //Initialization for all non-origin nodes
        for (int i = 0; i < number_of_nodes; ++i)
        {
            // not scanned
            m_node_status_array[i] = 0;
            m_node_label_cost[i] = _MAX_LABEL_COST;
            // pointer to previous NODE INDEX from the current label at current node and time
            m_link_predecessor[i] = -1;
            // pointer to previous NODE INDEX from the current label at current node and time
            m_node_predecessor[i] = -1;
            // comment out to speed up comuting
            ////m_label_time_array[i] = 0;
            ////m_label_distance_array[i] = 0;
        }

        std::map<int, int> node_id_visit_mapping;   // used to prevent negative cost loop in the path

        // int internal_debug_flag = 0;
        if (NodeForwardStarArray[origin_node].OutgoingLinkSize == 0)
            return 0;

        //Initialization for origin node at the preferred departure time, at departure time, cost = 0, otherwise, the delay at origin node
        m_label_time_array[origin_node] = 0;
        m_node_label_cost[origin_node] = 0.0;
        //Mark:	m_label_distance_array[origin_node] = 0.0;

        // Peiheng, 02/05/21, duplicate initialization, remove them later
        // pointer to previous NODE INDEX from the current label at current node and time
        m_link_predecessor[origin_node] = -1;
        // pointer to previous NODE INDEX from the current label at current node and time
        m_node_predecessor[origin_node] = -1;

        SEList_clear();
        SEList_push_back(origin_node);

        std::map<int, int>::iterator it, it_begin, it_end;
        int from_node, to_node;
        int link_sqe_no;
        double new_time = 0;
        double new_distance = 0;
        double new_to_node_cost = 0;
        int tempFront;
        while (!(m_ListFront == -1))   //SEList_empty()
        {
            // from_node = SEList_front();
            // SEList_pop_front();  // remove current node FromID from the SE list

            from_node = m_ListFront;//pop a node FromID for scanning

            node_id_visit_mapping[from_node] = 1;  // map visit

            tempFront = m_ListFront;
            m_ListFront = m_SENodeList[m_ListFront];
            m_SENodeList[tempFront] = -1;

            m_node_status_array[from_node] = 2;

            if (dtalog.log_path() >= 2 || local_debugging_flag)
            {
                dtalog.output() << "SP:scan SE node: " << g_node_vector[from_node].node_id << " with "
                    << NodeForwardStarArray[from_node].OutgoingLinkSize << " outgoing link(s). " << endl;
            }
            //scan all outbound nodes of the current node

            int pred_link_seq_no = m_link_predecessor[from_node];

            // for each link (i,j) belong A(i)
            for (int i = 0; i < NodeForwardStarArray[from_node].OutgoingLinkSize; ++i)
            {
                to_node = NodeForwardStarArray[from_node].OutgoingNodeNoArray[i];
                link_sqe_no = NodeForwardStarArray[from_node].OutgoingLinkNoArray[i];

                if (dtalog.log_path() >= 2 || local_debugging_flag)
                    dtalog.output() << "SP:  checking outgoing node " << g_node_vector[to_node].node_id << endl;

                // if(map (pred_link_seq_no, link_sqe_no) is prohibitted )
                //     then continue; //skip this is not an exact solution algorithm for movement

                if (g_node_vector[from_node].prohibited_movement_size >= 1)
                {
                    if (pred_link_seq_no >= 0)
                    {
                        string	movement_string;
                        string ib_link_id = g_link_vector[pred_link_seq_no].link_id;
                        string ob_link_id = g_link_vector[link_sqe_no].link_id;
                        movement_string = ib_link_id + "->" + ob_link_id;

                        if (g_node_vector[from_node].m_prohibited_movement_string_map.find(movement_string) != g_node_vector[from_node].m_prohibited_movement_string_map.end())
                        {
                            dtalog.output() << "prohibited movement " << movement_string << " will not be used " << endl;
                            continue;
                        }
                    }
                }

                //remark: the more complicated implementation can be found in paper Shortest Path Algorithms In Transportation Models: Classical and Innovative Aspects
                //	A note on least time path computation considering delays and prohibitions for intersection movements

                if (m_link_outgoing_connector_zone_seq_no_array[link_sqe_no] >= 0)
                {
                    if (m_link_outgoing_connector_zone_seq_no_array[link_sqe_no] != origin_zone)
                    {
                        // filter out for an outgoing connector with a centriod zone id different from the origin zone seq no
                        continue;
                    }
                }

                //very important: only origin zone can access the outbound connectors,
                //the other zones do not have access to the outbound connectors

                // Mark				new_time = m_label_time_array[from_node] + pLink->travel_time_per_period[tau];
                // Mark				new_distance = m_label_distance_array[from_node] + pLink->length;
                float additional_cost = 0;

                if (g_link_vector[link_sqe_no].RT_travel_time > 1)  // used in real time routing only
                {
                    additional_cost = g_link_vector[link_sqe_no].RT_travel_time;

                    //if (g_link_vector[link_sqe_no].RT_travel_time > 999)
                    //    continue; //skip this link due to closure
                }


                new_to_node_cost = m_node_label_cost[from_node] + m_link_genalized_cost_array[link_sqe_no] + additional_cost;

                if (dtalog.log_path() || local_debugging_flag)
                {
                    dtalog.output() << "SP:  checking from node " << g_node_vector[from_node].node_id
                        << "  to node " << g_node_vector[to_node].node_id << " cost = " << new_to_node_cost << endl;
                }

                bool b_visit_flag = false;

                if (node_id_visit_mapping.find(to_node) != node_id_visit_mapping.end())
                {
                    b_visit_flag = true;
                }


                if (b_visit_flag==false && new_to_node_cost < m_node_label_cost[to_node]) // we only compare cost at the downstream node ToID at the new arrival time t
                {

                    if (dtalog.log_path() || local_debugging_flag)
                    {
                        dtalog.output() << "SP:  updating node: " << g_node_vector[to_node].node_id << " current cost:" << m_node_label_cost[to_node]
                            << " new cost " << new_to_node_cost << endl;
                    }

                    // update cost label and node/time predecessor
                    // m_label_time_array[to_node] = new_time;
                    // m_label_distance_array[to_node] = new_distance;
                    m_node_label_cost[to_node] = new_to_node_cost;
                    // pointer to previous physical NODE INDEX from the current label at current node and time
                    m_node_predecessor[to_node] = from_node;
                    // pointer to previous physical NODE INDEX from the current label at current node and time
                    m_link_predecessor[to_node] = link_sqe_no;

                    if (dtalog.log_path() || local_debugging_flag)
                    {
                        dtalog.output() << "SP: add node " << g_node_vector[to_node].node_id << " new cost:" << new_to_node_cost
                            << " into SE List " << g_node_vector[to_node].node_id << endl;
                    }

                    // deque updating rule for m_node_status_array
                    if (m_node_status_array[to_node] == 0)
                    {
                        ///// SEList_push_back(to_node);
                        ///// begin of inline block
                        if (m_ListFront == -1)  // start from empty
                        {
                            m_ListFront = to_node;
                            m_ListTail = to_node;
                            m_SENodeList[to_node] = -1;
                        }
                        else
                        {
                            m_SENodeList[m_ListTail] = to_node;
                            m_SENodeList[to_node] = -1;
                            m_ListTail = to_node;
                        }
                        ///// end of inline block

                        m_node_status_array[to_node] = 1;
                    }

                    if (m_node_status_array[to_node] == 2)
                    {
                        /////SEList_push_front(to_node);
                        ///// begin of inline block
                        if (m_ListFront == -1)  // start from empty
                        {
                            m_SENodeList[to_node] = -1;
                            m_ListFront = to_node;
                            m_ListTail = to_node;
                        }
                        else
                        {
                            m_SENodeList[to_node] = m_ListFront;
                            m_ListFront = to_node;
                        }
                        ///// end of inline block

                        m_node_status_array[to_node] = 1;
                    }
                }
            }
        }

        if (dtalog.log_path() || local_debugging_flag)
        {
            dtalog.output() << "SPtree at iteration k = " << iteration_k << " origin node: "
                << g_node_vector[origin_node].node_id << endl;

            //Initialization for all non-origin nodes
            for (int i = 0; i < p_assignment->g_number_of_nodes; ++i)
            {
                int node_pred_id = -1;
                int node_pred_no = m_node_predecessor[i];

                if (node_pred_no >= 0)
                    node_pred_id = g_node_vector[node_pred_no].node_id;

                if (m_node_label_cost[i] < 9999)
                {
                    dtalog.output() << "SP node: " << g_node_vector[i].node_id << " label cost " << m_node_label_cost[i] << "time "
                        << m_label_time_array[i] << "node_pred_id " << node_pred_id << endl;
                }
            }
        }

        if (d_node_no >= 1)
            return m_node_label_cost[d_node_no];
        else
            return 0;  // one to all shortest pat
    }

    ////////// link based SE List

    int m_LinkBasedSEListFront;
    int m_LinkBasedSEListTail;
    int* m_LinkBasedSEList;  // dimension: number of links

    // SEList: Scan List implementation: the reason for not using STL-like template is to avoid overhead associated pointer allocation/deallocation
    void LinkBasedSEList_clear()
    {
        m_LinkBasedSEListFront = -1;
        m_LinkBasedSEListTail = -1;
    }

    void LinkBasedSEList_push_front(int link)
    {
        if (m_LinkBasedSEListFront == -1)  // start from empty
        {
            m_LinkBasedSEList[link] = -1;
            m_LinkBasedSEListFront = link;
            m_LinkBasedSEListTail = link;
        }
        else
        {
            m_LinkBasedSEList[link] = m_LinkBasedSEListFront;
            m_LinkBasedSEListFront = link;
        }
    }

    void LinkBasedSEList_push_back(int link)
    {
        if (m_LinkBasedSEListFront == -1)  // start from empty
        {
            m_LinkBasedSEListFront = link;
            m_LinkBasedSEListTail = link;
            m_LinkBasedSEList[link] = -1;
        }
        else
        {
            m_LinkBasedSEList[m_LinkBasedSEListTail] = link;
            m_LinkBasedSEList[link] = -1;
            m_LinkBasedSEListTail = link;
        }
    }

    bool LinkBasedSEList_empty()
    {
        return(m_LinkBasedSEListFront == -1);
    }

    int LinkBasedSEList_front()
    {
        return m_LinkBasedSEListFront;
    }

    void LinkBasedSEList_pop_front()
    {
        int tempFront = m_LinkBasedSEListFront;
        m_LinkBasedSEListFront = m_LinkBasedSEList[m_LinkBasedSEListFront];
        m_LinkBasedSEList[tempFront] = -1;
    }
};

