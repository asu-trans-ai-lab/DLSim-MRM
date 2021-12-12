
#ifndef GUARD_DTA_H
#define GUARD_DTA_H
#define BUILD_EXE //self-use

constexpr auto _MAX_LABEL_COST = 1.0e+15;
constexpr auto _INFO_ZONE_ID = 100000;

constexpr auto _MAX_AGNETTYPES = 7; //because of the od demand store format,the MAX_demandtype must >=g_DEMANDTYPES.size()+1;
constexpr auto _MAX_TIMEPERIODS = 1; // time period set to 4: mid night, morning peak, mid-day and afternoon peak;
constexpr auto _MAX_MEMORY_BLOCKS = 100;

constexpr auto _MAX_LINK_SIZE_IN_A_PATH = 5000;		// lu
constexpr auto _MAX_LINK_SIZE_FOR_A_NODE = 500;

constexpr auto _MAX_TIMESLOT_PerPeriod = 100; // max 96 15-min slots per day
constexpr auto _default_saturation_flow_rate = 1530;

constexpr auto MIN_PER_TIMESLOT = 15;
constexpr auto _simulation_discharge_period = 1;

/* make sure we change the following two parameters together*/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
constexpr auto number_of_seconds_per_interval = 0.25;  // consistent with the cell length of 7 meters
constexpr auto number_of_simu_interval_reaction_time = 4;  // reaction time as 1 second, 4 simu intervals, CAV: 0.5 seconds

constexpr auto number_of_simu_intervals_in_min = 240; // 60/0.25 number_of_seconds_per_interval

/* number_of_seconds_per_interval should satisify the ratio of 60/number_of_seconds_per_interval is an integer*/

// Linear congruential generator
constexpr auto LCG_a = 17364;
constexpr auto LCG_c = 0;
constexpr auto LCG_M = 65521;  // it should be 2^32, but we use a small 16-bit number to save memory

// FILE* g_pFileOutputLog = nullptr;
extern void g_OutputModelFiles(int mode);
template <typename T>
T* Allocate1DDynamicArray(int nRows)
{
    T* dynamicVector;

    dynamicVector = new (std::nothrow) T[nRows]();

    if (dynamicVector == NULL)
    {
        exit(1);

    }
    return dynamicVector;
}

template <typename T>
void Deallocate1DDynamicArray(T* dVector, int nRows)
{
    if (!dVector)
        return;
    delete[] dVector;
}


template <typename T>
T** Allocate2DDynamicArray(int nRows, int nCols)
{
    T** dynamicArray;

    dynamicArray = new (std::nothrow) T * [nRows];

    if (!dynamicArray)
    {
        dtalog.output() << "Error: insufficient memory.";
        g_ProgramStop();
    }

    for (int i = 0; i < nRows; ++i)
    {
        dynamicArray[i] = new (std::nothrow) T[nCols];

        if (!dynamicArray[i])
        {
            dtalog.output() << "Error: insufficient memory.";
            g_ProgramStop();
        }
    }

    return dynamicArray;
}

template <typename T>
void Deallocate2DDynamicArray(T** dArray, int nRows, int nCols)
{
    if (!dArray)
        return;

    for (int x = 0; x < nRows; ++x)
        delete[] dArray[x];

    delete[] dArray;
}

template <typename T>
T*** Allocate3DDynamicArray(int nX, int nY, int nZ)
{
    T*** dynamicArray = new (std::nothrow) T * *[nX];

    if (!dynamicArray)
    {
        dtalog.output() << "Error: insufficient memory.";
        g_ProgramStop();
    }

    for (int x = 0; x < nX; ++x)
    {
        if (x % 1000 == 0)
        {
            dtalog.output() << "allocating 3D memory for " << x << endl;
        }

        dynamicArray[x] = new (std::nothrow) T * [nY];

        if (!dynamicArray[x])
        {
            dtalog.output() << "Error: insufficient memory.";
            g_ProgramStop();
        }

        for (int y = 0; y < nY; ++y)
        {
            dynamicArray[x][y] = new (std::nothrow) T[nZ];
            if (!dynamicArray[x][y])
            {
                dtalog.output() << "Error: insufficient memory.";
                g_ProgramStop();
            }
        }
    }

    for (int x = 0; x < nX; ++x)
        for (int y = 0; y < nY; ++y)
            for (int z = 0; z < nZ; ++z)
                dynamicArray[x][y][z] = 0;

    return dynamicArray;
}

template <typename T>
void Deallocate3DDynamicArray(T*** dArray, int nX, int nY)
{
    if (!dArray)
        return;

    for (int x = 0; x < nX; ++x)
    {
        for (int y = 0; y < nY; ++y)
            delete[] dArray[x][y];

        delete[] dArray[x];
    }

    delete[] dArray;
}

template <typename T>
T**** Allocate4DDynamicArray(int nM, int nX, int nY, int nZ)
{
    T**** dynamicArray = new (std::nothrow) T * **[nX];

    if (!dynamicArray)
    {
        dtalog.output() << "Error: insufficient memory.";
        g_ProgramStop();
    }

    if (nM == 0 || nX == 0 || nY == 0 || nZ == 0)
    {
        dtalog.output() << "allocating 4D memory but size = 0 in 1 dimension.";
        g_ProgramStop();
    }
    for (int m = 0; m < nM; ++m)
    {
        if (m % 1000 == 0)
        {

            dtalog.output() << "allocating 4D memory for " << m << " zones,"
                << "nX=" << nX << ","
                << "nY=" << nY << ","
                << "nZ=" << nZ << endl;

        }

        dynamicArray[m] = new (std::nothrow) T * *[nX];

        if (!dynamicArray[m])
        {
            dtalog.output() << "Error: insufficient memory.";
            g_ProgramStop();
        }

        for (int x = 0; x < nX; ++x)
        {
            dynamicArray[m][x] = new (std::nothrow) T * [nY];

            if (!dynamicArray[m][x])
            {
                dtalog.output() << "Error: insufficient memory.";
                g_ProgramStop();
            }

            for (int y = 0; y < nY; ++y)
            {
                dynamicArray[m][x][y] = new (std::nothrow) T[nZ];
                if (!dynamicArray[m][x][y])
                {
                    dtalog.output() << "Error: insufficient memory.";
                    g_ProgramStop();
                }
            }
        }
    }

    return dynamicArray;
}

template <typename T>
void Deallocate4DDynamicArray(T**** dArray, int nM, int nX, int nY)
{
    if (!dArray)
        return;

    for (int m = 0; m < nM; ++m)
    {
        for (int x = 0; x < nX; ++x)
        {
            for (int y = 0; y < nY; ++y)
                delete[] dArray[m][x][y];

            delete[] dArray[m][x];
        }

        delete[] dArray[m];
    }

    delete[] dArray;
}


class CDemand_Period {
public:
    CDemand_Period() : demand_period{ 0 }, starting_time_slot_no{ 0 }, ending_time_slot_no{ 0 }, m_RandomSeed{ 101 }
    {
    }

    int get_time_horizon_in_min()
    {
        return (ending_time_slot_no - starting_time_slot_no) * 15;
    }

    unsigned int m_RandomSeed;

    float GetRandomRatio()
    {
        //m_RandomSeed is automatically updated.
        m_RandomSeed = (LCG_a * m_RandomSeed + LCG_c) % LCG_M;

        return float(m_RandomSeed) / LCG_M;
    }

    void compute_cumulative_profile(int starting_slot_no, int ending_slot_no)
    {
        float total_ratio = 0;
        for (int s = starting_slot_no; s < ending_slot_no; s++)
        {
            total_ratio += departure_time_ratio[s];
        }
        if (total_ratio < 0.000001)
            total_ratio = 0.000001;

        cumulative_departure_time_ratio[starting_slot_no] = 0;
        float cumulative_ratio = 0;
        for (int s = starting_slot_no; s < ending_slot_no; s++)
        {
            cumulative_ratio += departure_time_ratio[s] / total_ratio;
            cumulative_departure_time_ratio[s] = cumulative_ratio;
            dtalog.output() << "cumulative profile ratio at slot  " << s << " = " << cumulative_departure_time_ratio[s] << endl;
        }
        dtalog.output() << "final cumulative profile ratio" << cumulative_departure_time_ratio[ending_slot_no - 1] << endl;

    }
    int get_time_slot_no()
    {
        float r = GetRandomRatio();
        for (int s = starting_time_slot_no; s < ending_time_slot_no; s++)
        {
            if (r < cumulative_departure_time_ratio[s])
                return s;
        }
        return starting_time_slot_no;  // first time slot as the default value
    }

    float departure_time_ratio[_MAX_TIMESLOT_PerPeriod];
    float cumulative_departure_time_ratio[_MAX_TIMESLOT_PerPeriod];

    string demand_period;
    int starting_time_slot_no;
    int ending_time_slot_no;
    string time_period;
    int demand_period_id;
};

class CAgent_type {
public:
    CAgent_type() : agent_type_no{ 1 }, value_of_time{ 1 }, time_headway_in_sec{ 1 }, real_time_information{ 0 }, access_speed{ 2 }, access_distance_lb{ 0.0001 }, access_distance_ub{ 4 }, acecss_link_k{ 4 }
    {
    }

    int agent_type_no;
    // dollar per hour
    float value_of_time;
    // link type, product consumption equivalent used, for travel time calculation
    float PCE;
    float time_headway_in_sec;
    int real_time_information;
    string agent_type;
    string display_code;

    string access_node_type;
    float access_speed;

    float access_distance_lb;
    float access_distance_ub;
    int acecss_link_k;


};

class CLinkType
{
public:
    CLinkType() : link_type{ 1 }, number_of_links{ 0 }, traffic_flow_code{ 0 }
    {
    }


    int link_type;
    int number_of_links;
    int traffic_flow_code;

    string link_type_name;
    string type_code;
};

class CColumnPath {
public:
    CColumnPath() : path_node_vector{ nullptr }, path_link_vector{ nullptr }, path_seq_no{ 0 },
        path_switch_volume{ 0 }, path_volume{ 0 }, path_travel_time{ 0 }, path_distance{ 0 }, path_toll{ 0 },
        path_gradient_cost{ 0 }, path_gradient_cost_difference{ 0 }, path_gradient_cost_relative_difference{ 0 }, subarea_output_flag{ 1 }, measurement_flag{ 0 }
    {
    }

    ~CColumnPath()
    {
        if (m_node_size >= 1)
        {
            delete[] path_node_vector;
            delete[] path_link_vector;
        }
    }

    // Peiheng, 02/02/21, consider using move later
    void AllocateVector(int node_size, const int* node_vector, int link_size, const int* link_vector, bool backwardflag = false)
    {
        m_node_size = node_size;
        m_link_size = link_size;
        // dynamic array
        path_node_vector = new int[node_size];
        path_link_vector = new int[link_size];

        if (backwardflag)
        {
            // copy backward
            for (int i = 0; i < m_node_size; ++i)
                path_node_vector[i] = node_vector[m_node_size - 1 - i];

            for (int i = 0; i < m_link_size; ++i)
                path_link_vector[i] = link_vector[m_link_size - 1 - i];
        }
        else
        {
            // copy forward
            for (int i = 0; i < m_node_size; ++i)
                path_node_vector[i] = node_vector[i];

            for (int i = 0; i < m_link_size; ++i)
                path_link_vector[i] = link_vector[i];
        }
    }

    // Peiheng, 02/02/21, consider using move later
    void AllocateVector(const std::vector<int>& node_vector, const std::vector<int>& link_vector, bool backwardflag = false)
    {
        m_node_size = node_vector.size();
        m_link_size = link_vector.size();
        // dynamic array
        path_node_vector = new int[m_node_size];
        path_link_vector = new int[m_link_size];

        if (backwardflag)
        {
            // copy backward
            for (int i = 0; i < m_node_size; ++i)
                path_node_vector[i] = node_vector[m_node_size - 1 - i];

            for (int i = 0; i < m_link_size; ++i)
                path_link_vector[i] = link_vector[m_link_size - 1 - i];
        }
        else
        {
            // copy forward
            for (int i = 0; i < m_node_size; ++i)
                path_node_vector[i] = node_vector[i];

            for (int i = 0; i < m_link_size; ++i)
                path_link_vector[i] = link_vector[i];
        }
    }

    int* path_node_vector;
    int* path_link_vector;

    std::vector<int> path_link_STL_vector;
    int path_seq_no;
    // path volume
    double path_volume;
    int subarea_output_flag;
    int measurement_flag;
    double path_switch_volume;
    double path_travel_time;
    double path_distance;
    double path_toll;
    // first order graident cost.
    double path_gradient_cost;
    // first order graident cost - least gradient cost
    double path_gradient_cost_difference;
    // first order graident cost - least gradient cost
    double path_gradient_cost_relative_difference;

    int m_node_size;
    int m_link_size;
    std::map <int, bool> diverted_vehicle_map;


    std::vector<int> agent_simu_id_vector;
};

class CAgentPath {
public:
    CAgentPath() : path_id{ 0 }, node_sum{ -1 }, travel_time{ 0 }, distance{ 0 }, volume{ 0 }
    {
    }

    int path_id;
    int node_sum;
    float travel_time;
    float distance;
    float volume;
    std::vector <int> path_link_sequence;


};

class CColumnVector {

public:
    // this is colletion of unique paths
    CColumnVector() : cost{ 0 }, time{ 0 }, distance{ 0 }, od_volume{ 0 }, bfixed_route{ false }, m_passing_sensor_flag{ -1 }, information_type{ 0 }, activity_agent_type_no{ 0 }
    {
    }

    float cost;
    float time;
    float distance;
    // od volume
    double od_volume;
    bool bfixed_route;
    int information_type;
    std::vector<int> activity_zone_no_vector;
    int activity_agent_type_no;



    int m_passing_sensor_flag;
    // first key is the sum of node id;. e.g. node 1, 3, 2, sum of those node ids is 6, 1, 4, 2 then node sum is 7.
    // Peiheng, 02/02/21, potential memory leak, fix it
    std::map <int, CColumnPath> path_node_sequence_map;
};

class CAgent_Column {
public:
    CAgent_Column() : cost{ 0 }
    {
    }

    float cost;
    float volume;
    float travel_time;
    float distance;

    int agent_id;
    int o_zone_id;
    int d_zone_id;
    int o_node_id;
    int d_node_id;

    string agent_type;
    string demand_period;

    vector<int> path_node_vector;
    vector<int> path_link_vector;
    vector<float> path_time_vector;
};

// event structure in this "event-based" traffic simulation

class DTAVehListPerTimeInterval
{
public:
    std::vector<int> m_AgentIDVector;
};



class CAgent_Simu
{
public:
    CAgent_Simu() : agent_vector_seq_no{ -1 }, path_toll{ 0 }, departure_time_in_min{ 0 }, m_bGenereated{ false }, m_bCompleteTrip{ false },
        path_travel_time_in_min{ 0 }, path_distance{ 0 }, diversion_flag{ 0 }, time_headway{ number_of_simu_interval_reaction_time }, PCE_unit_size{ 1 }
    {
    }

    ~CAgent_Simu()
    {

    }

    void AllocateMemory()
    {

        m_current_link_seq_no = 0;

        m_Veh_LinkArrivalTime_in_simu_interval.reserve(path_link_seq_no_vector.size());
        m_Veh_LinkDepartureTime_in_simu_interval.reserve(path_link_seq_no_vector.size());

        for (int i = 0; i < path_link_seq_no_vector.size(); ++i)
        {
            m_Veh_LinkArrivalTime_in_simu_interval.push_back(-1);
            m_Veh_LinkDepartureTime_in_simu_interval.push_back(-1);
        }

        m_path_link_seq_no_vector_size = path_link_seq_no_vector.size();
        departure_time_in_simu_interval = int(departure_time_in_min * 60.0 / number_of_seconds_per_interval + 0.5);  // round off

    }


    float GetRandomRatio()
    {
        // Peiheng, 02/02/21, m_RandomSeed is uninitialized
        //m_RandomSeed is automatically updated.
        m_RandomSeed = (LCG_a * m_RandomSeed + LCG_c) % LCG_M;

        return float(m_RandomSeed) / LCG_M;
    }

    int agent_vector_seq_no;
    float path_toll;
    float departure_time_in_min;
    int diversion_flag;
    bool m_bGenereated;
    bool m_bCompleteTrip;

    int agent_service_type;
    int demand_type;
    int agent_id;

    // for column pool index
    int at;
    int tau;
    int dest;

    int m_current_link_seq_no;
    int m_path_link_seq_no_vector_size;

    int departure_time_in_simu_interval;
    float arrival_time_in_min;
    float path_travel_time_in_min;
    float path_distance;

    unsigned int m_RandomSeed;

    // external input
    std::vector<int> path_link_seq_no_vector;
    std::vector<int>  m_Veh_LinkArrivalTime_in_simu_interval;
    std::vector<int>  m_Veh_LinkDepartureTime_in_simu_interval;
    int time_headway;  // in terms of simulation interval 
    int PCE_unit_size;  // the number of units: 
};



class Assignment {
public:
    // default is UE
    Assignment() : assignment_mode{ 0 }, g_number_of_memory_blocks{ 8 }, g_number_of_threads{ 1 }, g_link_type_file_loaded{ true }, g_agent_type_file_loaded{ false },
        total_demand_volume{ 0.0 }, g_column_pool{ nullptr }, g_number_of_in_memory_simulation_intervals{ 500 },
        g_number_of_column_generation_iterations{ 20 }, g_number_of_demand_periods{ 24 }, g_number_of_links{ 0 }, g_number_of_timing_arcs{ 0 },
        g_number_of_nodes{ 0 }, g_number_of_zones{ 0 }, g_number_of_agent_types{ 0 }, debug_detail_flag{ 1 }, path_output{ 1 }, trajectory_output{ 1 }, major_path_volume_threshold{ 0.000001 }, trajectory_sampling_rate{ 1.0 }, dynamic_link_performance_sampling_interval_in_min{ 60 }, dynamic_link_performance_sampling_interval_hd_in_min{ 15 }, trajectory_diversion_only{ 0 }, m_GridResolution{ 0.01 }
    {
    }

    ~Assignment()
    {
        if (g_column_pool)
            Deallocate4DDynamicArray(g_column_pool, g_number_of_zones, g_number_of_zones, g_number_of_agent_types);


        DeallocateLinkMemory4Simulation();
    }

    void InitializeDemandMatrix(int number_of_zones, int number_of_agent_types, int number_of_time_periods)
    {
        total_demand_volume = 0.0;
        g_number_of_zones = number_of_zones;
        g_number_of_agent_types = number_of_agent_types;

        g_column_pool = Allocate4DDynamicArray<CColumnVector>(number_of_zones, number_of_zones, max(1, number_of_agent_types), number_of_time_periods);

        for (int i = 0; i < number_of_zones; ++i)
        {
            g_origin_demand_array[i] = 0.0;
        }

        for (int i = 0; i < number_of_agent_types; ++i)
        {
            for (int tau = 0; tau < g_number_of_demand_periods; ++tau)
                total_demand[i][tau] = 0.0;
        }

        g_DemandGlobalMultiplier = 1.0f;
    }

    int get_in_memory_time(int t)
    {
        return t % g_number_of_in_memory_simulation_intervals;
    }

    void STTrafficSimulation();
    void STMesoTrafficSimulation();

    //OD demand estimation estimation
    void Demand_ODME(int OD_updating_iterations);
    void AllocateLinkMemory4Simulation();
    void UpdateRTPath(CAgent_Simu* pAgent);
    bool RTSP_RealTimeShortestPathFinding(int time_slot_no, int simu_interval_t);
    void DeallocateLinkMemory4Simulation();

    double m_GridResolution;
    int assignment_mode;
    int g_number_of_memory_blocks;
    int g_number_of_threads;
    int path_output;
    int trajectory_output;
    float trajectory_sampling_rate;
    int trajectory_diversion_only;
    int dynamic_link_performance_sampling_interval_in_min;
    float dynamic_link_performance_sampling_interval_hd_in_min;

    float major_path_volume_threshold;

    bool g_link_type_file_loaded;
    bool g_agent_type_file_loaded;

    float total_demand_volume;
    std::map<int, float> g_origin_demand_array;
    CColumnVector**** g_column_pool;

    // the data horizon in the memory
    int g_number_of_in_memory_simulation_intervals;
    int g_number_of_column_generation_iterations;
    int g_number_of_demand_periods;


    std::map<_int64, int> cell_id_mapping;  // this is used to mark if this cell_id has been identified or not;

    int g_number_of_links;
    int g_number_of_timing_arcs;
    int g_number_of_nodes;
    int g_number_of_zones;
    int g_number_of_agent_types;

    std::map<int, int> node_seq_no_2_info_zone_id_mapping;  // this is used to mark if this zone_id has been identified or not
    std::map<int, int> zone_seq_no_2_info_mapping;  // this is used to mark if this zone_id has been identified or not
    std::map<int, int> zone_seq_no_2_activity_mapping;  // this is used to mark if this zone_id has been identified or not

    int debug_detail_flag;

    // hash table, map external node number to internal node sequence no.
    std::map<int, int> g_node_id_to_seq_no_map;
    std::map<string, int> g_mvmt_key_to_link_no_map;
    // from integer to integer map zone_id to zone_seq_no
    std::map<int, int> g_zoneid_to_zone_seq_no_mapping;
    std::map<string, int> g_link_id_map;

    std::map<int, double> zone_id_X_mapping;
    std::map<int, double> zone_id_Y_mapping;

    std::vector<CDemand_Period> g_DemandPeriodVector;
    int g_LoadingStartTimeInMin;
    int g_LoadingEndTimeInMin;

    std::vector<CAgent_type> g_AgentTypeVector;
    std::map<int, CLinkType> g_LinkTypeMap;

    std::map<string, int> demand_period_to_seqno_mapping;
    std::map<string, int> agent_type_2_seqno_mapping;

    float total_demand[_MAX_AGNETTYPES][_MAX_TIMEPERIODS];
    float g_DemandGlobalMultiplier;

    // used in ST Simulation
    float** m_LinkOutFlowCapacity;  // per second interval for simplicity
    int** m_LinkOutFlowState;  // per second interval for simplicity


    // in min
    float** m_LinkTDWaitingTime;
    std::vector<float> m_LinkTotalWaitingTimeVector;;
    // number of simulation time intervals

    float** m_LinkCumulativeArrivalVector;
    float** m_LinkCumulativeDepartureVector;

    float* m_LinkCACount;  // CA, assign this value to m_LinkCumulativeArrivalVector at a given time in min
    float* m_LinkCDCount;  // CD

    int g_start_simu_interval_no;
    int g_number_of_simulation_intervals;
    // is shorter than g_number_of_simulation_intervals
    int g_number_of_loading_intervals_in_sec;
    // the data horizon in the memory in min
    int g_number_of_intervals_in_min;

    int g_number_of_intervals_in_sec;

};

extern Assignment assignment;

class CVDF_Period
{
public:
    CVDF_Period() : m{ 0.5 }, VOC{ 0 }, gamma{ 3.47f }, mu{ 1000 }, PHF{ 3 },
        alpha{ 0.15f }, beta{ 4 }, rho{ 1 }, preload{ 0 }, penalty{ 0 }, LR_price{ 0 }, LR_RT_price{ 0 }, marginal_base{ 1 },
        starting_time_slot_no{ 0 }, ending_time_slot_no{ 0 },
        cycle_length{ -1 }, red_time{ 0 }, effective_green_time{ 0 }, t0{ 0 }, t3{ 0 }, start_green_time{ -1 }, end_green_time{ -1 }
    {
        for (int at = 0; at < _MAX_AGNETTYPES; at++)
        {
            toll[at] = 0;
            pce[at] = 0;
        }

        for (int t = 0; t < _MAX_TIMESLOT_PerPeriod; ++t)
        {
            Queue[t] = 0;
            waiting_time[t] = 0;
            arrival_rate[t] = 0;
            discharge_rate[t] = 0;
            travel_time[t] = 0;
        }
    }

    ~CVDF_Period()
    {
    }

    float get_waiting_time(int relative_time_slot_no)
    {
        if (relative_time_slot_no >= 0 && relative_time_slot_no < _MAX_TIMESLOT_PerPeriod)
            return waiting_time[relative_time_slot_no];
        else
            return 0;
    }

    float PerformSignalVDF(float hourly_per_lane_volume, float red, float cycle_length)
    {
        float lambda = hourly_per_lane_volume;
        float mu = _default_saturation_flow_rate; //default saturation flow ratesa
        float s_bar = 1.0 / 60.0 * red * red / (2 * cycle_length); // 60.0 is used to convert sec to min
        float uniform_delay = s_bar / max(1 - lambda / mu, 0.1f);

        return uniform_delay;
    }

    double PerformBPR(double volume)
    {
        // take nonnegative values
        volume = max(0.0, volume);

        // Peiheng, 02/02/21, useless block
        if (volume > 1.0)
        {
            int debug = 1;
        }

        VOC = volume / max(0.00001, capacity);
        avg_travel_time = FFTT + FFTT * alpha * pow((volume + preload) / max(0.00001, capacity), beta);
        total_travel_time = (volume + preload) * avg_travel_time;
        marginal_base = FFTT * alpha * beta * pow((volume + preload) / max(0.00001, capacity), beta - 1);

        if (cycle_length > 1)
        {
            if (volume > 10)
                avg_travel_time += cycle_length / 60.0 / 2.0; // add additional random delay due to signalized intersection by 
            else
            {
                int debug = 1;
            }
        }
        return avg_travel_time;
        // volume --> avg_traveltime
    }

    // input period based volume
    double PerformBPR_X(double volume)
    {
        bValidQueueData = false;
        congestion_period_P = 0;

        float FFTT_in_hour = FFTT / 60.0;
        // convert avg_travel_time from unit of min to hour
        float avg_travel_time_in_hour = avg_travel_time / 60.0;

        // Step 1: Initialization
        int L = ending_time_slot_no - starting_time_slot_no;  // in 15 min slot
        if (L >= _MAX_TIMESLOT_PerPeriod - 1)
            return 0;

        for (int t = starting_time_slot_no; t <= ending_time_slot_no; ++t)
        {
            Queue[t] = 0;
            waiting_time[t] = 0;
            arrival_rate[t] = 0;
            discharge_rate[t] = mu / 2.0;
            travel_time[t] = FFTT_in_hour;
        }

        // avg_travel time should be per min
        // avg_travel_time = (FFTT_in_hour + avg_waiting_time)*60.0;

        //int L = ending_time_slot_no - starting_time_slot_no;  // in 15 min slot
        float mid_time_slot_no = starting_time_slot_no + L / 2.0;  // t1;  // we can discuss thi
        // Case 1: fully uncongested region
        if (volume <= L * mu / 2)
        {
            // still keep 0 waiting time for all time period
            congestion_period_P = 0;
        }
        else
        {
            // partially congested region
            //if (volume > L * mu / 2 ) // Case 2
            float P = 0;
            // unit: hour  // volume / PHF  is D, D/mu = congestion period
            congestion_period_P = volume / PHF / mu;
            P = congestion_period_P * 4; //unit: 15 time slot

            t0 = max(0.0, mid_time_slot_no - P / 2.0);
            t3 = min((double)_MAX_TIMESLOT_PerPeriod - 1, mid_time_slot_no + P / 2.0);

            // we need to recalculate gamma coefficient based on D/C ratio. based on assumption of beta = 4;
            // average waiting time based on eq. (32)
            // https://www.researchgate.net/publication/348488643_Introduction_to_connection_between_point_queue_model_and_BPR
            // w= gamma*power(D/mu,4)/(80*mu)
            // gamma = w*80*mu/power(D/mu,4)
            // avg_travel_time has been calculated based on standard BPR function, thus, the condition is that, PerformBPR() should be called before PerformBPR_X

            float uncongested_travel_time_in_hour = FFTT_in_hour;
            float w = avg_travel_time_in_hour - uncongested_travel_time_in_hour; //unit: hour
            // mu is read from the external link.csv file
            float D_over_mu = congestion_period_P;

            gamma = w * 120 * mu / pow(D_over_mu, 4);  // m =1/2, gamma = w80... when m = 3/4

            // derive t0 and t3 based on congestion duration p
            int t2 = m * (t3 - t0) + t0;
            //tt_relative is relative time
            for (int tt_relative = 0; tt_relative <= L; ++tt_relative)
            {
                //absolute time index
                int time_abs = starting_time_slot_no + tt_relative;
                if (time_abs < t0)
                {
                    //first uncongested phase with mu/2 as the approximate flow rates
                    waiting_time[time_abs] = 0;  // per hour
                    arrival_rate[time_abs] = mu / 2;
                    discharge_rate[time_abs] = mu / 2.0;
                    travel_time[time_abs] = FFTT_in_hour;  // per hour
                }

                if (time_abs >= t0 && time_abs <= t3 && t3 > t0)
                {
                    float t_ph = time_abs / (60.0 / MIN_PER_TIMESLOT);
                    float t0_ph = t0 / (60.0 / MIN_PER_TIMESLOT);
                    float t2_ph = t2 / (60.0 / MIN_PER_TIMESLOT);
                    float t3_ph = t3 / (60.0 / MIN_PER_TIMESLOT);

                    //second congested phase based on the gamma calculated from the dynamic eq. (32)
                    Queue[time_abs] = 1 / (4.0) * gamma * (t_ph - t0_ph) * (t_ph - t0_ph) * (t_ph - t3_ph) * (t_ph - t3_ph);
                    // unit is hour
                    waiting_time[time_abs] = 1 / (4.0 * mu) * gamma * (t_ph - t0_ph) * (t_ph - t0_ph) * (t_ph - t3_ph) * (t_ph - t3_ph);
                    arrival_rate[time_abs] = gamma * (t_ph - t0_ph) * (t_ph - t2_ph) * (t_ph - t3_ph) + mu;
                    discharge_rate[time_abs] = mu;
                    travel_time[time_abs] = FFTT_in_hour + waiting_time[time_abs]; // per hour
                    bValidQueueData = true;
                }

                if (time_abs > t3)
                {
                    //third uncongested phase with mu/2 as the approximate flow rates
                    waiting_time[time_abs] = 0;
                    arrival_rate[time_abs] = mu / 2;
                    discharge_rate[time_abs] = mu / 2.0;
                    travel_time[time_abs] = FFTT_in_hour;
                }
                // avg_waiting_time = gamma / (120 * mu)*pow(P, 4.0) *60.0;// avg_waiting_time  should be per min
                // dtalog() << avg_waiting_time << endl;
                // avg_travel_time = (FFTT_in_hour + avg_waiting_time)*60.0; // avg_travel time should be per min
            }
        }

        return avg_travel_time;
    }

    double m;
    // we should also pass uncongested_travel_time as length/(speed_at_capacity)
    double VOC;
    //updated BPR-X parameters
    double gamma;
    double mu;
    //peak hour factor
    double PHF;
    //standard BPR parameter
    double alpha;
    double beta;
    double preload;
    double toll[_MAX_AGNETTYPES];
    double pce[_MAX_AGNETTYPES];
    double penalty;
    double LR_price[_MAX_AGNETTYPES];
    double LR_RT_price[_MAX_AGNETTYPES];;
    string allowed_uses;


    double rho;
    double marginal_base;
    // in 15 min slot
    int starting_time_slot_no;
    int ending_time_slot_no;
    float cycle_length;
    float red_time;
    float effective_green_time;
    int start_green_time;
    int end_green_time;

    int t0, t3;

    bool bValidQueueData;
    string period;

    double capacity;
    double FFTT;

    double congestion_period_P;
    // inpput
    double volume;

    //output
    double avg_delay;
    double avg_travel_time = 0;
    double avg_waiting_time = 0;
    double total_travel_time = 0;

    // t starting from starting_time_slot_no if we map back to the 24 hour horizon
    float Queue[_MAX_TIMESLOT_PerPeriod];
    float waiting_time[_MAX_TIMESLOT_PerPeriod];
    float arrival_rate[_MAX_TIMESLOT_PerPeriod];

    float discharge_rate[_MAX_TIMESLOT_PerPeriod];
    float travel_time[_MAX_TIMESLOT_PerPeriod];
};

class CLink
{
public:
    // construction
    CLink() :main_node_id{ -1 }, obs_count{ -1 }, upper_bound_flag{ 0 }, est_count_dev{ 0 }, free_speed{ 0 },
        BWTT_in_simulation_interval{ 100 }, zone_seq_no_for_outgoing_connector{ -1 }, number_of_lanes{ 1 }, lane_capacity{ 1999 },
        length{ 1 }, free_flow_travel_time_in_min{ 1 }, link_spatial_capacity{ 100 }, 
        timing_arc_flag{ false }, traffic_flow_code{ 0 }, spatial_capacity_in_vehicles{ 999999 }, link_type{ 2 }, subarea_id{ -1 }, RT_flow_volume{ 0 },
        cell_type{ -1 }, saturation_flow_rate{ 1800 }, dynamic_link_reduction_start_time_slot_no{ 99999 }, b_automated_generated_flag {false}
    {
        for (int tau = 0; tau < _MAX_TIMEPERIODS; ++tau)
        {
            flow_volume_per_period[tau] = 0;
            queue_length_perslot[tau] = 0;
            travel_time_per_period[tau] = 0;
            TDBaseTT[tau] = 0;
            TDBaseCap[tau] = 0;
            TDBaseFlow[tau] = 0;
            TDBaseQueue[tau] = 0;
            //cost_perhour[tau] = 0;

            for (int at = 0; at < _MAX_AGNETTYPES; ++at)
                volume_per_period_per_at[tau][at] = 0;
        }

        //for (int tau = 0; tau < _MAX_TIMESLOT_PerPeriod; ++tau)
        //{
        //    RT_travel_time_vector[tau] = -1;
        //    RT_speed_vector[tau] = -1;

        //}
    }

    ~CLink()
    {
    }

    // Peiheng, 02/05/21, useless block
    void free_memory()
    {
    }

    void Calculatedynamic_VDFunction();
    void Calculatedynamic_RTVDFunction();


    float get_VOC_ratio(int tau)
    {
        return (flow_volume_per_period[tau] + TDBaseFlow[tau]) / max(0.00001, TDBaseCap[tau]);
    }

    float get_speed(int tau)
    {
        return length / max(travel_time_per_period[tau], 0.0001) * 60;  // per hour
    }

    void calculate_marginal_cost_for_agent_type(int tau, int agent_type_no, float PCE_agent_type)
    {
        // volume * dervative
        // BPR_term: volume * FFTT * alpha * (beta) * power(v/c, beta-1),

//        travel_marginal_cost_per_period[tau][agent_type_no] = VDF_period[tau].marginal_base * PCE_agent_type;
    }

    float get_generalized_first_order_gradient_cost_of_second_order_loss_for_agent_type(int tau, int agent_type_no)
    {
        // *60 as 60 min per hour
        float generalized_cost = travel_time_per_period[tau] + VDF_period[tau].penalty + VDF_period[tau].toll[agent_type_no] / assignment.g_AgentTypeVector[agent_type_no].value_of_time * 60;

        // system optimal mode or exterior panalty mode
        //if (assignment.assignment_mode == 4)
        //    generalized_cost += travel_marginal_cost_per_period[tau][agent_type_no];

        return generalized_cost;
    }

    int main_node_id;

    float obs_count;
    int upper_bound_flag;
    float est_count_dev;

    int BWTT_in_simulation_interval;
    int zone_seq_no_for_outgoing_connector;

    int number_of_lanes;
    double lane_capacity;
    double saturation_flow_rate;

    std::map <int, int> m_link_pedefined_capacity_map;  // per sec
// please remove this redundent data structure 
    float dynamic_link_capacity[_MAX_TIMESLOT_PerPeriod];
    int dynamic_link_reduction_start_time_slot_no;
    std::map <int, bool> dynamic_link_closure_map;
    std::map <int, string> dynamic_link_closure_type_map;

    double length;
    double free_flow_travel_time_in_min;
    double free_speed;

    double cost;
    double link_spatial_capacity;

    bool timing_arc_flag;
    int traffic_flow_code;
    int spatial_capacity_in_vehicles;
    int time_to_be_released;

    // 1. based on BPR.

    int link_seq_no;

    std::map<int, int> capacity_reduction_map;
    
    string link_id;
    string geometry;
    bool AllowAgentType(string agent_type, int tau)
    {
        if (VDF_period[tau].allowed_uses.size() == 0 || VDF_period[tau].allowed_uses == "all")  // if the allowed_uses is empty then all types are allowed.
            return true;
        else
        {
            if (VDF_period[tau].allowed_uses.find(agent_type) != string::npos)  // otherwise, only an agent type is listed in this "allowed_uses", then this agent type is allowed to travel on this link
                return true;
            else
            {
                return false;
            }


        }
    }

    int from_node_seq_no;
    int to_node_seq_no;
    int link_type;
    bool b_automated_generated_flag;

    int cell_type;
    string mvmt_txt_id;
    string path_code_str;
    string tmc_corridor_name;
    string link_type_name;
    string link_type_code;

    float PCE;
    float fftt;

    CVDF_Period VDF_period[_MAX_TIMEPERIODS];

    double TDBaseTT[_MAX_TIMEPERIODS];
    double TDBaseCap[_MAX_TIMEPERIODS];
    double TDBaseFlow[_MAX_TIMEPERIODS];
    double TDBaseQueue[_MAX_TIMEPERIODS];

    int type;

    //static
    //float flow_volume;
    //float travel_time;

    int subarea_id;
    double flow_volume_per_period[_MAX_TIMEPERIODS];
    double RT_flow_volume;
    double background_flow_volume_per_period[_MAX_TIMEPERIODS];

    double  volume_per_period_per_at[_MAX_TIMEPERIODS][_MAX_AGNETTYPES];

    double  queue_length_perslot[_MAX_TIMEPERIODS];  // # of vehicles in the vertical point queue
    double travel_time_per_period[_MAX_TIMEPERIODS];
    double RT_travel_time;

    std::map<int, float> RT_travel_time_vector;
    std::map<int, float> RT_speed_vector;
    //	double  travel_marginal_cost_per_period[_MAX_TIMEPERIODS][_MAX_AGNETTYPES];

    int number_of_periods;

    //std::vector <SLinkMOE> m_LinkMOEAry;
    //beginning of simulation data

    //toll related link
    //int m_TollSize;
    //Toll *pTollVector;  // not using SLT here to avoid issues with OpenMP

    // for discrete event simulation
    // link-in queue of each link
    std::list<int> EntranceQueue;
    // link-out queue of each link
    std::list<int> ExitQueue;
    int current_driving_AgentID;
    int win_count;
    int lose_count;

};




class CNode
{
public:
    CNode() : zone_id{ -1 }, zone_org_id{ -1 }, prohibited_movement_size{ 0 }, node_seq_no{ -1 }, subarea_id{ -1 }, is_activity_node{ 0 }, is_information_zone{ 0 }, agent_type_no{ -1 }
    {
    }

    //int accessible_node_count;

    int zone_id;
    __int64 cell_id;
    string cell_str;
    // original zone id for non-centriod nodes
    int zone_org_id;

    string node_type;
    string agent_type_str;
    int subarea_id;
    int prohibited_movement_size;
    // sequence number
    int node_seq_no;

    //external node number
    int node_id;

    int is_activity_node;
    int is_information_zone;
    int agent_type_no;

    double x;
    double y;

    std::vector<int> m_outgoing_link_seq_no_vector;
    std::vector<int> m_incoming_link_seq_no_vector;

    std::vector<int> m_to_node_seq_no_vector;
    std::map<int, int> m_to_node_2_link_seq_no_map;

    std::map<string, int> m_prohibited_movement_string_map;
    // for simulation
    std::map<int, int> next_link_for_resource_request;
};

class CInfoCell {
public:
    __int64 cell_id;
    string cell_str;
    std::vector<GDPoint> m_ShapePoints;

    void CreateCell(double x, double y, double grid_resolution)
    {

        __int64 xi;
        xi = floor(x / grid_resolution);

        __int64 yi;
        yi = floor(y / grid_resolution);

        double left, right, top, bottom;

        left = xi * grid_resolution;
        right = (xi + 1) * grid_resolution;

        top = (yi + 1) * grid_resolution;
        bottom = (yi)*grid_resolution;

        GDPoint	pt0, pt1, pt2, pt3, pt4;

        pt0.x = left; 	pt0.y = top;
        pt1.x = right; 	pt1.y = top;
        pt2.x = right; 	pt2.y = bottom;
        pt3.x = left; 	pt3.y = bottom;
        pt4.x = left; 	pt4.y = top;


        m_ShapePoints.push_back(pt0);
        m_ShapePoints.push_back(pt1);
        m_ShapePoints.push_back(pt2);
        m_ShapePoints.push_back(pt3);
        m_ShapePoints.push_back(pt4);

    }

};

class CTMC_Corridor_Info {
public:
    CTMC_Corridor_Info()
    {
    }

    void reset()
    {
        total_VMT = 0;
        total_VHT = 0;
        total_VDT = 0;
        lowest_speed = 9999;
        highest_speed = -1;
        link_count = 0;
    }

    double get_avg_speed()
    {
        return total_VMT / max(0.001, total_VHT);  //miles per hour
    }
    double total_VMT;
    double total_VHT;
    double total_VDT;

    double avg_speed;
    double lowest_speed;
    double highest_speed;
    double total_congestion_duration;
    int link_count;

    std::map<int, int> road_sequence_map;

};



class COZone
{
public:
    COZone() : obs_production{ -1 }, obs_attraction{ -1 },
        est_production{ -1 }, est_attraction{ -1 },
        est_production_dev{ 0 }, est_attraction_dev{ -1 }, total_demand{ 0 }, b_real_time_information{ false }
    {
    }

    bool b_real_time_information;

    float obs_production;
    float obs_attraction;

    float est_production;
    float est_attraction;

    float est_production_dev;
    float est_attraction_dev;

    // 0, 1,
    int zone_seq_no;
    // external zone id // this is origin zone
    int zone_id;
    int node_seq_no;
    float total_demand;
    float obs_production_upper_bound_flag;
    float obs_attraction_upper_bound_flag;

    std::vector<int> m_activity_node_vector;
};

extern std::vector<COZone> g_zone_vector;

class CAGBMAgent
{
public:
    int agent_id;
    int income;
    int gender;
    int vehicle;
    int purpose;
    int flexibility;
    float preferred_arrival_time;
    float travel_time_in_min;
    float free_flow_travel_time;
    int from_zone_seq_no;
    int to_zone_seq_no;
    int type;
    int time_period;
    int k_path;
    float volume;
    float arrival_time_in_min;
};



extern std::vector<CNode> g_node_vector;
extern std::vector<CLink> g_link_vector;

extern std::map<int, DTAVehListPerTimeInterval> g_AgentTDListMap;
extern vector<CAgent_Simu*> g_agent_simu_vector;;

extern std::map<string, CTMC_Corridor_Info> g_tmc_corridor_vector;
extern std::map<string, CInfoCell> g_info_cell_map;

#endif
