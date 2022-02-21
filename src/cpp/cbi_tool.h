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

class CTMC_Link
{
public:

    CTMC_Link()
    {
        b_with_sensor_speed_data = false;
        for (int t = 0; t < MAX_TIMEINTERVAL_PerDay; ++t)
        {
            speed_sum[t] = 0;
            avg_speed[t] = -1;
            speed_lowest[t] = 99999;

            volume_sum[t] = 0;
            avg_volume[t] = -1;

            b_volume_data_available_flag[t] = false;
            speed_count[t] = 0;
            volume_count[t] = 0;
        }

    }

    ~CTMC_Link()
    {
    }

    float speed_lowest[MAX_TIMEINTERVAL_PerDay];
    float volume_per_hour_per_lane[MAX_TIMEINTERVAL_PerDay];
    float speed_count[MAX_TIMEINTERVAL_PerDay];
    float volume_count[MAX_TIMEINTERVAL_PerDay];



    float get_highest_speed()
    {
        float highest_speed = 0;

        for (int t_in_min = 6 * 60; t_in_min < 20 * 60; t_in_min += 5)
        {

            float avg_speed = record_avg_speed(t_in_min);

            if (avg_speed > highest_speed)
                highest_speed = avg_speed;
        }

        return highest_speed;
    }

    float scan_highest_speed_and_vc(float& FD_vcutoff, float free_speed, float & highest_speed)
    {
        highest_speed = 0;
        for (int t_in_min = 6 * 60; t_in_min < 20 * 60; t_in_min += 5)
        {

            float avg_speed = record_avg_speed(t_in_min);

            if (avg_speed > highest_speed)
                highest_speed = avg_speed;
        }

        return FD_vcutoff / max(1.0f, free_speed) * highest_speed;  // return the final cutoff speed
    }

    float check_feasible_range(float input_value, float default_value, float lower_bound, float upper_bound)
    {
        float output_value = input_value;

        if (input_value < lower_bound || input_value >upper_bound)
        {
            output_value = default_value;
        }
        return output_value;
    }


    float scan_congestion_duration(int peak_no, float starting_time_in_hour, float ending_time_in_hour, float assign_period_t2_peak_in_hour, float& FD_vcutoff, CLink* p_link,  /* input*/
        float& obs_t0_in_hour, float& obs_t3_in_hour, float& obs_P_in_hour,
        float& V, float& peak_hour_volume, float&D, float& VOC_ratio, float& DOC_ratio,
        float& mean_speed_BPR, float &mean_speed_QVDF, float& highest_speed,  float& t2_speed,
        float& plf, float &qdf, float &Q_n, float& Q_s, float& Q_cd, float& Q_cp)
    {

       Q_n = 1.24;
       Q_s = 4;
       Q_cd = 1;
       Q_cp = 0.24;

       
        obs_t0_in_hour = -1;
        obs_t3_in_hour = -1;
        obs_P_in_hour = 0;
        V = 0;

        D = 0;
        DOC_ratio = 0;
        VOC_ratio = 0;
        mean_speed_BPR = 0;
        mean_speed_QVDF = 0;
        highest_speed = 0;

        int obs_t0_in_interval = -1;
        int obs_t3_in_interval = -1;
        int obs_ts_in_interval = starting_time_in_hour * 12;
        int obs_te_in_interval = ending_time_in_hour * 12;


        float L = ending_time_in_hour - starting_time_in_hour;
        plf = 1.0/max(0.1f, L);
        qdf = 1.0/max(0.1f, L);

        float total_speed_value = 0;
        float total_speed_count = 0;

        // step 1: record higest speed
        for (int t_in_min = 6 * 60; t_in_min < 20 * 60; t_in_min += 5)
        {

            float avg_speed = record_avg_speed(t_in_min);

            if (avg_speed > highest_speed)
                highest_speed = avg_speed;
        }

        int t_mid = assign_period_t2_peak_in_hour * 12;

        V = 0;

        int t_lowest_speed = -1;
        double lowest_speed = 99999;
        t2_speed = avg_speed[t_mid];

        // step 3: compute volume V
        if (avg_speed[t_mid] > 0)
        {
            for (int t_in_min = starting_time_in_hour * 60; t_in_min < ending_time_in_hour * 60; t_in_min += 5)
            {

                float avg_speed = record_avg_speed(t_in_min);
                float volume = get_avg_volume(t_in_min,p_link,avg_speed, highest_speed);
                V += volume / 12; // 12 5-min interval per hour

                if (avg_speed < lowest_speed)
                {
                    t_lowest_speed = t_in_min / 5;
                    lowest_speed = avg_speed;
                }

            }

        }

        if (avg_speed[t_mid] < 0)
            return 0;

        // below is congested period analysis

        bool b_t0_found_flag = false;
        bool b_t3_found_flag = false;
        //step 5: 
        if (avg_speed[t_mid] < FD_vcutoff)  // exit if the speed is higher than v_congestion_cutoff
        {
            obs_t3_in_interval = t_mid;
            obs_t3_in_hour = t_mid * 1.0 / 12.0;
            b_t3_found_flag = true;

            obs_t0_in_interval = t_mid;
            b_t0_found_flag = true;
            obs_t0_in_hour = t_mid * 1.0 / 12.0;

        }

        for (int t = t_mid + 1; t <= obs_te_in_interval; t += 1)  // move forward from mid t
        {
            if (avg_speed[t] < 1)
            {
                int i_no_data = 1;
                //                dtalog.output() << " Error: TMC_link " << p_link->tmc_code.c_str() << " has no data at timt t " << t << g_time_coding(t * 5).c_str() << endl;;
                break;
            }

           if (avg_speed[t] > FD_vcutoff)  // exit if the speed is higher than v_congestion_cutoff
            {
                break;
            }

           obs_t3_in_interval = t;
           obs_t3_in_hour = t * 1.0 / 12.0;
           b_t3_found_flag = true;


        }


        for (int t = t_mid - 1; t >= obs_ts_in_interval; t -= 1) // move backward from t_mid 
        {
            if (avg_speed[t] < 1)
            {
                int i_no_data = 1;
                //dtalog.output() << " Error: TMC link " << p_link->tmc_code.c_str() << " has no data at timt t =" << t * 5 << " min" << endl;;
                break;
            }
            if (avg_speed[t] > FD_vcutoff)
            {
                break;
            }

            obs_t0_in_interval = t;
            b_t0_found_flag = true;
            obs_t0_in_hour = t * 1.0 / 12.0;

        }

        // consider peak hour
        int peak_hour_t0_in_interval = t_mid - 60 / 2 / 5;  // - 30 min
        int peak_hour_t3_in_interval = t_mid + 60 / 2 / 5;  //+ 30 min

        total_speed_count = 0;  // initial values
        total_speed_value = 0;


        obs_P_in_hour = 0;  // no compromise 
        D = 0;
        DOC_ratio = 0;

        peak_hour_volume = 0;

        if (p_link->link_id == "201102BA")
            int debug_flag = 1;

        // consider voluem in peak hour
        for (int t = peak_hour_t0_in_interval; t < peak_hour_t3_in_interval; t += 1) // move between congestion duration per interval
        {
            if (avg_speed[t] >= 1)
            {
                total_speed_value += avg_speed[t];
                float volume = get_avg_volume(t, p_link, avg_speed[t], highest_speed);

                if (volume < 0)
                {
                    int idebug = 1;
                    get_avg_volume(t, p_link, avg_speed[t], highest_speed);
                }
                peak_hour_volume += volume / 12; // 12 5-min interval per hour
                total_speed_count++;
            }

        }

        mean_speed_BPR = total_speed_value / max(1.0f, total_speed_count);

        plf = peak_hour_volume / max(1.0f, V);
        VOC_ratio = peak_hour_volume / max(1.0, p_link->lane_capacity);  // unit: demand: # of vehicles, lane_capacity # of vehicles per hours: dc ratio has a unit of hour, but it is different from P

        mean_speed_QVDF = FD_vcutoff;  

        if (b_t0_found_flag == false && b_t3_found_flag == true)
        {
            obs_t0_in_hour = assign_period_t2_peak_in_hour;
            b_t0_found_flag = true; // reset
        }

        if (b_t0_found_flag == true && b_t3_found_flag == false)
        {
            obs_t3_in_hour = assign_period_t2_peak_in_hour;
            b_t3_found_flag = true;
        }

        if (b_t0_found_flag == false || b_t3_found_flag == false)  // two boundaries are not found or P < 1;
        {
            //uncongested states

            // consider peak hour
            obs_t0_in_interval = t_mid - 60 / 2 / 5;  // - 30 min
            obs_t3_in_interval = t_mid + 60 / 2 / 5;  //+ 30 min


            obs_P_in_hour = 0;  // no compromise 
            D = 0;
            DOC_ratio = 0;

            qdf = 1.0 / L;

        }else
        {  // congested states
            D = 0;
            total_speed_count = 0;  // initial values
            total_speed_value = 0;
            double lowest_speed = FD_vcutoff;
            int time_int_with_lowest_speed = obs_t0_in_interval;

            for (int t = obs_t0_in_interval; t <= obs_t3_in_interval; t += 1) // move between congestion duration per interval
            {
                if (avg_speed[t] >= 1)
                {
                    total_speed_value += avg_speed[t];
                    float volume = get_avg_volume(t, p_link, avg_speed[t], highest_speed);
                    D += volume / 12; // 12 5-min interval per hour
                    total_speed_count++;

                    if (avg_speed[t] < lowest_speed)
                    {
                        lowest_speed = avg_speed[t];
                        time_int_with_lowest_speed = t;
                    }
                }

            }

            // test
            obs_P_in_hour = (obs_t3_in_hour - obs_t0_in_hour);  // congestion duration P

                DOC_ratio = D / max(1.0, p_link->lane_capacity);  // unit: demand: # of vehicles, lane_capacity # of vehicles per hours: dc ratio has a unit of hour, but it is different from P
                VOC_ratio = max(VOC_ratio, DOC_ratio);

            mean_speed_QVDF = total_speed_value / max(1.0f, total_speed_count);

            t2_speed = lowest_speed;  // if we use a pure second order model, we should consider t2= 2/3(t3-t0)+ t0

            if (obs_P_in_hour > 2.5)
            {
                int idebug = 1;
            }
            // calibration
            double exact_qdf = D / max(1.0f, V);
            double exact_plf = peak_hour_volume / max(1.0f, V);
            double exact_bound = 1.0 / L;
            plf = qdf = max(exact_bound, max( exact_plf, exact_qdf));  // pure qdf

            //revise qdf output as the result of plf 
            // P = cd * (D/C) ^n  --> log (P) = log(cd) + n *Log (D/C)
            double part1 = (log(obs_P_in_hour) - log(Q_cd));
            double part2 = log(DOC_ratio);

            if (abs(part2) < 0.000001)
                part2 = 0.00001;

            Q_n = part1 / part2;  // assume Q_cd is fixed at 1 o other values close to 1
            
            if (Q_n < 1.001)
            {
                Q_n = 1.124; // default, to ensure the mu is decreasing pattern as a function of D/C
                //Cd = P / (D / C) ^ n
                // Peiheng, 02/22/22, g++ will treat pow(DOC_ratio, Q_n) as double but not clang++
                Q_cd = obs_P_in_hour / max(0.0001, (double)pow(DOC_ratio, Q_n));
            }

            //vc / vt2 - 1 = cp * (P)^s, --> cp = [ vc/vt2 - 1] / (P^s)  // assume s is fixed
            // Peiheng, 02/22/22, g++ will treat pow(obs_P_in_hour, Q_s) as double but not clang++
            Q_cp = (FD_vcutoff / max(0.0001f, t2_speed) - 1.0) / max(0.00001, (double)pow(obs_P_in_hour, Q_s));
            //backward derivation
        }

        // modified with Mohammad A. 02/15/2022 
        Q_n = check_feasible_range(Q_n, 1.24, 1,1.5);
        Q_s = check_feasible_range(Q_s, 1, 0.5, 4);
        Q_cd = check_feasible_range(Q_cd, 1, 0.5, 2);
        Q_cp = check_feasible_range(Q_cp, 0.2, 0.0, 2);
        double largest_DOC_ratio = max (1.0f, L);
        DOC_ratio = check_feasible_range(DOC_ratio, 0.5, 0.0, largest_DOC_ratio);
        plf = check_feasible_range(plf, 1.0 / L, 0.0, 1);
        qdf = check_feasible_range(qdf, 1.0 / L, 0.0, 1);

        obs_P_in_hour = check_feasible_range(obs_P_in_hour, 0.0, 0.0, 10);
        t2_speed = check_feasible_range(t2_speed, FD_vcutoff, 0.0, 200);
        return obs_P_in_hour;
    }

    string tmc_code;

    int matching_link_no;
    bool b_with_sensor_speed_data;

    // construction
    void add_speed_sensor_data(int day_no, int time_in_min, float speed_value, float volume_pl_value)
    {
        b_with_sensor_speed_data = true;
        int t = time_in_min / 5;

        if (t >= 0 && t < MAX_TIMEINTERVAL_PerDay)
        {
            speed_sum[t] += speed_value;
            speed_count[t] += 1;

            if (speed_value < speed_lowest[t])
                speed_lowest[t] = speed_value;

        }

        if (t >= 0 && t < MAX_TIMEINTERVAL_PerDay && volume_pl_value >=0)
        {
            b_volume_data_available_flag[t] = true;
            volume_sum[t] += volume_pl_value;
            volume_count[t] += 1;
        }

    }

    float record_avg_speed(int time_in_min)
    {
        int t = time_in_min / 5;

        if (speed_count[t] == 0)
            return -1;

        if (t >= 0 && t < MAX_TIMEINTERVAL_PerDay)
        {

            avg_speed[t] = speed_sum[t] / max(1.0f, speed_count[t]);
            return avg_speed[t];
        }

        return -1;
    }



    float get_avg_speed(int time_in_min)
    {
        int t = time_in_min / 5;
        if (t >= 0 && t < MAX_TIMEINTERVAL_PerDay)
        {

            return speed_sum[t] / max(1.0f, speed_count[t]);

        }
        return -1;
    }



    float get_avg_volume(int time_in_min, CLink* p_link, float avg_speed, float highest_speed)
    {
        int t = time_in_min / 5;
            
        if (volume_count[t] >0 && t >= 0 && t < MAX_TIMEINTERVAL_PerDay)
        {

            return volume_sum[t] / max(1.0f, volume_count[t]) *12;  // 5 min to hourly volume

        }
        else  // default 
        {
            return p_link->get_volume_from_speed(avg_speed, highest_speed, p_link->lane_capacity);
        }

    }
    float get_avg_speed_15min(int time_in_min)
    {
        int t = time_in_min / 5;

        float total_speed_value = 0;
        int total_speed_count = 0;

        for (int tt = 0; tt < 3; tt++)
        {

            if (t + tt >= 0 && t + tt < MAX_TIMEINTERVAL_PerDay)
            {
                total_speed_value += speed_sum[t + tt] / max(1.0f, speed_count[t + tt]);
                total_speed_count++;
            }
        }

        return total_speed_value / max(1, total_speed_count);

    }

    float get_avg_hourly_speed(int time_in_min)
    {
        int t = time_in_min / 5;

        float total_speed_value = 0;
        int total_speed_count = 0;

        for (int tt = 0; tt < 12; tt++)
        {

            if (t + tt >= 0 && t + tt < MAX_TIMEINTERVAL_PerDay)
            {
                total_speed_value += speed_sum[t + tt] / max(1.0f, speed_count[t + tt]);
                total_speed_count++;
            }
        }

        return total_speed_value / max(1, total_speed_count);
    }



    float get_avg_hourly_volume(int time_in_min)
    {
        int t = time_in_min / 5;

        float total_volume_value = 0;
        int total_volume_count = 0;

        for (int tt = 0; tt < 12; tt++)
        {

            if (t + tt >= 0 && t + tt < MAX_TIMEINTERVAL_PerDay)
            {
                total_volume_value += volume_sum[t + tt] / max(1.0f, volume_count[t + tt]);
                total_volume_count++;
            }
        }

        return total_volume_value;
    }

    float get_lowest_speed(int time_in_min)
    {
        int t = time_in_min / 5;
        if (t >= 0 && t < MAX_TIMEINTERVAL_PerDay)
        {

            return speed_lowest[t];

        }

    }


    float reference_speed;

    float avg_speed[MAX_TIMEINTERVAL_PerDay];
    float speed_sum[MAX_TIMEINTERVAL_PerDay];
    float avg_volume[MAX_TIMEINTERVAL_PerDay];
    float volume_sum[MAX_TIMEINTERVAL_PerDay];
    bool   b_volume_data_available_flag[MAX_TIMEINTERVAL_PerDay];

};

extern std::vector<CTMC_Link> g_TMC_vector;
std::map<int, int> g_dayDataMap;
int g_dayofweek(int y, int m, int d)
{
    static int t[] = { 0, 3, 2, 5, 0, 3,
                       5, 1, 4, 6, 2, 4 };
    y -= m < 3;
    return (y + y / 4 - y / 100 +
        y / 400 + t[m - 1] + d) % 7;
}

int g_dayofyear(int y, int m, int d)
{
    return (m - 1) * 31 + d - 1;

}
float g_measurement_tstamp_parser(string str, int& day_of_week_flag, int& day_of_year)
{

    int string_lenghth = str.length();

    //ASSERT(string_lenghth < 100);

    const char* string_line = str.data(); //string to char*

    int char_link_distance_VDF = strlen(string_line);

    char ch, buf_hh[32] = { 0 }, buf_mm[32] = { 0 }, buf_ss[32] = { 0 };

    int yyyy, month, day;
    char  hh1, hh2, mm1, mm2;
    float hhf1, hhf2, mmf1, mmf2;
    float global_minute = 0;
    float dd = 0, hh = 0, mm = 0;
    int i = 11;  //2019-01-01 06:00:00,
    int buffer_i = 0, buffer_k = 0, buffer_j = 0;
    int num_of_colons = 0;
    int num_of_underscore = 0;

    int char_link_distance_VDF_yymmdd = 10;
    i = 0;

    //    sscanf(string_line, "%d/%d/%d", &month, &day, &yyyy );
    sscanf(string_line, "%d-%d-%d", &yyyy, &month, &day);



    day_of_week_flag = g_dayofweek(yyyy, month, day);
    day_of_year = g_dayofyear(yyyy, month, day);

    g_dayDataMap[day_of_year] = month * 100 + day;


    /// 
    i = 11;

    while (i < char_link_distance_VDF)
    {
        ch = string_line[i++];

        if (num_of_colons == 0 && ch != '_' && ch != ':') //input to buf_ddhhmm until we meet the colon
        {
            buf_hh[buffer_i++] = ch;
        }
        else if (num_of_colons == 1 && ch != ':') //start the Second "SS"
        {
            buf_mm[buffer_k++] = ch;
        }
        else if (num_of_colons == 2 && ch != ':') //start the Millisecond "sss"
        {
            buf_ss[buffer_j++] = ch;
        }

        if (i == char_link_distance_VDF) //start a new time string
        {
            //HHMM, 0123
            hh1 = buf_hh[0]; //read each first
            hh2 = buf_hh[1];
            mm1 = buf_mm[0];
            mm2 = buf_mm[1];

            hhf1 = ((float)hh1 - 48); //convert a char to a float
            hhf2 = ((float)hh2 - 48);
            mmf1 = ((float)mm1 - 48);
            mmf2 = ((float)mm2 - 48);

            dd = 0;
            hh = hhf1 * 10 * 60 + hhf2 * 60;
            mm = mmf1 * 10 + mmf2;

            global_minute = dd + hh + mm;

            //initialize the parameters
            buffer_i = 0;
            buffer_k = 0;
            buffer_j = 0;
            num_of_colons = 0;
            break;
        }

        if (ch == ':')
            num_of_colons += 1;
    }

    return global_minute;
}
bool Assignment::map_tmc_reading()
{
    // step 1: read measurement.csv
    CCSVParser parser_measurement;


    if (parser_measurement.OpenCSVFile("Reading.csv", true))
    {
        int count = 0;
        while (parser_measurement.ReadRecord())  // if this line contains [] mark, then we will also read field headers.
        {
            string tmc;

            parser_measurement.GetValueByFieldName("tmc_code", tmc);
            string measurement_tstamp;

            float tmc_reference_speed = 0;
            float speed = -1;
            float volume_pl = -1;
            bool bMatched_flag = false;
            parser_measurement.GetValueByFieldName("measurement_tstamp", measurement_tstamp, false);

            float global_time;
            int day_of_week_flag = 0;
            int day_of_year = 0;

            if (measurement_tstamp.size() < 18)
                continue; // skip empty lines 

            global_time = g_measurement_tstamp_parser(measurement_tstamp, day_of_week_flag, day_of_year);

            if (day_of_week_flag == 0 && day_of_week_flag == 6)
                continue;

            parser_measurement.GetValueByFieldName("speed", speed);
            parser_measurement.GetValueByFieldName("volume_pl", volume_pl, false);

            
            parser_measurement.GetValueByFieldName("reference_speed", tmc_reference_speed, false);

            string ROADNAME;
            parser_measurement.GetValueByFieldName("ROADNAME", ROADNAME, false);

            if (assignment.m_TMClink_map.find(tmc) == assignment.m_TMClink_map.end())
            {
                CTMC_Link tmc_link;
                tmc_link.tmc_code = tmc;
                assignment.m_TMClink_map[tmc] = g_TMC_vector.size();
                g_TMC_vector.push_back(tmc_link);

            }

            int index = assignment.m_TMClink_map[tmc];

            g_TMC_vector[index].add_speed_sensor_data(day_of_year, global_time, speed, volume_pl);


            if (count % 100000 == 0)
            {
                dtalog.output() << "reading " << count / 100000 << "00k TMC data items" << endl;

            }
            count++;
        }
        parser_measurement.CloseCSVFile();

        dtalog.output() << "reading data for " << g_TMC_vector.size() << " TMC links." << endl;
        return true;
    }
    else
    {
        dtalog.output() << "reading data with " << g_TMC_vector.size() << " TMC links." << endl;
        g_program_stop();
    }

    return false;

}
void g_output_tmc_file()
{
    // if 
    FILE* p_file_tmc_link = fopen("link_cbi_summary.csv", "w");


    if (p_file_tmc_link != NULL)
    {
        fprintf(p_file_tmc_link, "link_id,tmc,tmc_corridor_name,tmc_corridor_id,tmc_road_order,tmc_road_sequence,tmc_road,tmc_direction,tmc_intersection,tmc_highest_speed,link_no,from_node_id,to_node_id,link_type,");
        fprintf(p_file_tmc_link, "link_type_code,FT,AT,vdf_code,nlanes,link_distance_VDF,free_speed,capacity,k_critical,vcutoff,highest_speed,vcutoff_updated,vcutoff_ratio,v_critical_s3,");

        for (int tau = 0; tau < min((size_t)3, assignment.g_DemandPeriodVector.size()); tau++)
        {

            fprintf(p_file_tmc_link, "%s_t0,%s_t3,%s_V,%s_peak_hour_volume,%s_D,%s_VC_ratio,%s_DC_ratio,%s_P,%s_vc/vt2-1,%s_vf_delay_index,%s_vc_delay_index,%s_speed_ph,%s_queue_speed,%s_vt2,%s_qdf,%s_Q_n,%s_Q_cp,%s_Q_alpha,%s_Q_beta,%s_mV,%s_mD,%s_mDC_ratio,%s_mP,%s_mv_QVDF,%s_mvt2_QVDF,%s_m_mu_QVDF,%s_m_gamma_QVDF,%s_m_peak_hour_volume,%s_mVC_ratio,%s_mv_BPR,",
               assignment.g_DemandPeriodVector[tau].demand_period.c_str(),
               assignment.g_DemandPeriodVector[tau].demand_period.c_str(),
               assignment.g_DemandPeriodVector[tau].demand_period.c_str(),
               assignment.g_DemandPeriodVector[tau].demand_period.c_str(),
               assignment.g_DemandPeriodVector[tau].demand_period.c_str(),

               assignment.g_DemandPeriodVector[tau].demand_period.c_str(),
               assignment.g_DemandPeriodVector[tau].demand_period.c_str(),
               assignment.g_DemandPeriodVector[tau].demand_period.c_str(),
               assignment.g_DemandPeriodVector[tau].demand_period.c_str(),
               assignment.g_DemandPeriodVector[tau].demand_period.c_str(),
               
               assignment.g_DemandPeriodVector[tau].demand_period.c_str(),
               assignment.g_DemandPeriodVector[tau].demand_period.c_str(),
               assignment.g_DemandPeriodVector[tau].demand_period.c_str(),
                assignment.g_DemandPeriodVector[tau].demand_period.c_str(),
                assignment.g_DemandPeriodVector[tau].demand_period.c_str(),
               
                assignment.g_DemandPeriodVector[tau].demand_period.c_str(),
               assignment.g_DemandPeriodVector[tau].demand_period.c_str(),
               assignment.g_DemandPeriodVector[tau].demand_period.c_str(),
                assignment.g_DemandPeriodVector[tau].demand_period.c_str(),
                assignment.g_DemandPeriodVector[tau].demand_period.c_str(),
                
                assignment.g_DemandPeriodVector[tau].demand_period.c_str(),
                assignment.g_DemandPeriodVector[tau].demand_period.c_str(),
               assignment.g_DemandPeriodVector[tau].demand_period.c_str(),
               assignment.g_DemandPeriodVector[tau].demand_period.c_str(),
               assignment.g_DemandPeriodVector[tau].demand_period.c_str(), 
               assignment.g_DemandPeriodVector[tau].demand_period.c_str(),
               
                assignment.g_DemandPeriodVector[tau].demand_period.c_str(),
                assignment.g_DemandPeriodVector[tau].demand_period.c_str(), 
                assignment.g_DemandPeriodVector[tau].demand_period.c_str(), 
                assignment.g_DemandPeriodVector[tau].demand_period.c_str()
            );

        }


        fprintf(p_file_tmc_link, "geometry, tmc_geometry,");

        for (int t = 6 * 60; t < 20 * 60; t += 60)
        {
            int hour = t / 60;
            int minute = t - hour * 60;
            fprintf(p_file_tmc_link, "vh%02d,", hour);
        }

        for (int t = 6 * 60; t < 20 * 60; t += 60)
        {
            int hour = t / 60;
            int minute = t - hour * 60;
            fprintf(p_file_tmc_link, "m_vh%02d,", hour);
        }

        fprintf(p_file_tmc_link, "evhMAE,evhMAPE,evhRMSE,");

        for (int t = 6 * 60; t < 20 * 60; t += 60)
        {
            int hour = t / 60;
            int minute = t - hour * 60;
            fprintf(p_file_tmc_link, "qh%02d,", hour);
        }

        for (int t = 6 * 60; t < 20 * 60; t += 60)
        {
            int hour = t / 60;
            int minute = t - hour * 60;
            fprintf(p_file_tmc_link, "vhr%02d,", hour);
        }

        for (int t = 6 * 60; t < 20 * 60; t += 5)
        {
            int hour = t / 60;
            int minute = t - hour * 60;

            fprintf(p_file_tmc_link, "v%02d:%02d,", hour, minute);
        }

        for (int t = 6 * 60; t < 20 * 60; t += 15)
        {
            int hour = t / 60;
            int minute = t - hour * 60;

            fprintf(p_file_tmc_link, "v%02d:%02d,", hour, minute);
        }
        for (int t = 6 * 60; t < 20 * 60; t += 5)
        {
            int hour = t / 60;
            int minute = t - hour * 60;

            fprintf(p_file_tmc_link, "q%02d:%02d,", hour, minute);
        }



        fprintf(p_file_tmc_link, "\n");

        std::map<__int64, int> TMC_long_id_mapping;  // this is used to mark if this cell_id has been identified or not

        //// sort data records
        for (int i = 0; i < g_link_vector.size(); i++)
        {


            if (g_link_vector[i].tmc_code.size() > 0)
            {

                __int64 TMC_long_key = (g_link_vector[i].tmc_corridor_id * 10000 + g_link_vector[i].tmc_road_sequence) * 10 + g_link_vector[i].link_seq_no;
                TMC_long_id_mapping[TMC_long_key] = g_link_vector[i].link_seq_no;
            }
        }


        std::map<__int64, int>::iterator it;

        for (it = TMC_long_id_mapping.begin(); it != TMC_long_id_mapping.end(); ++it)
        {
            int i = it->second;
            float highest_speed = g_link_vector[i].free_speed;

            if (assignment.m_TMClink_map.find(g_link_vector[i].tmc_code) != assignment.m_TMClink_map.end())
            {
                int tmc_index = assignment.m_TMClink_map[g_link_vector[i].tmc_code];

                if (g_TMC_vector[tmc_index].b_with_sensor_speed_data == false)
                    continue;
                
                highest_speed = g_TMC_vector[tmc_index].get_highest_speed();
            }
            else
            {
                continue; // no tmc record data
            }


            float free_speed = g_link_vector[i].free_speed;

            if (g_link_vector[i].lane_capacity > 5000)  // skip connectors
                continue;

            g_link_vector[i].update_kc(g_link_vector[i].free_speed);


            fprintf(p_file_tmc_link, "%s,%s,%s,%d,%d,%d,%s,%s,%s,%.2f,%d,%d,%d,%d,%s,%d,%d,%s,%.1f,%f,%f,%f,%f,%f,",
                g_link_vector[i].link_id.c_str(),
                g_link_vector[i].tmc_code.c_str(),
                g_link_vector[i].tmc_corridor_name.c_str(),
                g_link_vector[i].tmc_corridor_id,
                g_link_vector[i].tmc_road_order,
                g_link_vector[i].tmc_road_sequence,
                g_link_vector[i].tmc_road.c_str(),
                g_link_vector[i].tmc_direction.c_str(),
                g_link_vector[i].tmc_intersection.c_str(),
                highest_speed,
                g_link_vector[i].link_seq_no,
                g_node_vector[g_link_vector[i].from_node_seq_no].node_id,
                g_node_vector[g_link_vector[i].to_node_seq_no].node_id,
                g_link_vector[i].link_type,
                g_link_vector[i].link_type_code.c_str(),
                g_link_vector[i].FT,
                g_link_vector[i].AT,
                g_link_vector[i].vdf_code.c_str(),
                g_link_vector[i].number_of_lanes,
                g_link_vector[i].link_distance_VDF,
                g_link_vector[i].free_speed,
                g_link_vector[i].lane_capacity,
                g_link_vector[i].k_critical,
                g_link_vector[i].v_congestion_cutoff

            );


            if (assignment.m_TMClink_map.find(g_link_vector[i].tmc_code) == assignment.m_TMClink_map.end())
            {
                fprintf(p_file_tmc_link, "\n");
                continue;

            }

            if (assignment.m_TMClink_map.find(g_link_vector[i].tmc_code) != assignment.m_TMClink_map.end())
            {
                int tmc_index = assignment.m_TMClink_map[g_link_vector[i].tmc_code];


                float highest_speed = 0;
                CLink* p_link = &(g_link_vector[i]);


                //if (p_link->link_id == "38644")
                //{
                //    int idebug = 1;
                //}
                float updated_vc = g_TMC_vector[tmc_index].scan_highest_speed_and_vc(p_link->v_congestion_cutoff, p_link->free_speed, highest_speed);
                p_link->v_congestion_cutoff = updated_vc;
                p_link->update_kc(free_speed);
                fprintf(p_file_tmc_link, "%f,%f,%f,", highest_speed, p_link->v_congestion_cutoff, p_link->v_congestion_cutoff / max(1.0f, highest_speed));
                fprintf(p_file_tmc_link, "%f,", p_link->v_critical);

                int analysis_hour_flag[25];

                for (int hour = 0; hour <= 24; hour++)
                    analysis_hour_flag[hour] = 0;


                //p_link->VDF_period[0].queue_demand_factor = 0.46;  //AM
                //p_link->VDF_period[1].queue_demand_factor = 0.24;  //MD
                //p_link->VDF_period[2].queue_demand_factor = 0.50;  //PM

                for (int tau = 0; tau < min ((size_t)3, assignment.g_DemandPeriodVector.size()); tau++)
                {

                    float assign_period_start_time_in_hour = assignment.g_DemandPeriodVector[tau].starting_time_slot_no * MIN_PER_TIMESLOT / 60.0;
                    float assign_period_end_time_in_hour = assignment.g_DemandPeriodVector[tau].ending_time_slot_no * MIN_PER_TIMESLOT / 60.0;
                    float assign_period_t2_peak_in_hour = assignment.g_DemandPeriodVector[tau].t2_peak_in_hour;

                    for (int hour = assign_period_start_time_in_hour; hour <= assign_period_end_time_in_hour; hour++)
                        analysis_hour_flag[hour] = 1;



                    //float scan_congestion_duration(int peak_no, float starting_time_in_hour, float ending_time_in_hour
                     // , float& FD_vcutoff, CLink * p_link,  /* input*/
                     //    float& obs_t0_in_hour, float& obs_t3_in_hour, float& obs_P,
                     //    float& V, float & peak_hour_volume, float& D, float& VOC_ratio,float& DOC_ratio,
                     //    float& mean_speed_BPR, float& mean_speed_QVDF, float& highest_speed, float& t2_speed, 
                     //float & plf, float &qdf )


                    float obs_t0_in_hour, obs_t3_in_hour, obs_P;
                    
                    obs_t0_in_hour = obs_t3_in_hour = obs_P = 0;
                    
                    float V, D, VOC_ratio, DOC_ratio;
                    V = D = VOC_ratio = DOC_ratio = 0;
                    
                    float peak_hour_volume = 0;
                    float mean_speed_BPR = 0;
                    float mean_speed_QVDF = 0;
                    float t2_speed = 0;
                    float plf= 1;
                    float qdf = 1;

                    float Q_n = 1;
                    float Q_s = 1;
                    float Q_cd = 1;
                    float Q_cp = 1;


                    obs_P = g_TMC_vector[tmc_index].scan_congestion_duration(tau,
                        assign_period_start_time_in_hour,
                        assign_period_end_time_in_hour,
                        assign_period_t2_peak_in_hour,
                        p_link->v_congestion_cutoff, p_link,
                        obs_t0_in_hour, obs_t3_in_hour, obs_P,
                        V, peak_hour_volume, D, VOC_ratio, DOC_ratio,
                        mean_speed_BPR, mean_speed_QVDF, highest_speed, t2_speed,
                        plf, qdf,
                        Q_n, Q_s, Q_cd, Q_cp);

                    qdf = plf;

                    // replacing the model data  based on observations only if the VDF setting has not been setup yet
                    // thus, if we have assignment results in a post processing step, we will use the assignment volume directly 
                        // location based parameters 
                    double Q_alpha = 8.0 / 15 * Q_cp * pow(Q_cd, Q_s);
                    double Q_beta = Q_n * Q_s;

                    Q_alpha = g_TMC_vector[tmc_index].check_feasible_range(Q_alpha, 0.27, 0.01, 1);
                    Q_alpha = g_TMC_vector[tmc_index].check_feasible_range(Q_beta, 1.14, 0.5, 5);

                    if (p_link->link_id == "201065AB")
                    {
                        int idebug = 1;

                    }
                    if (p_link->VDF_period[tau].queue_demand_factor < 0)
                    {

                        p_link->VDF_period[tau].sa_volume = V * p_link->number_of_lanes;  // make the data ready for sensitivity analysis
                        p_link->VDF_period[tau].peak_load_factor = plf;
                        p_link->VDF_period[tau].queue_demand_factor = qdf;
                        p_link->VDF_period[tau].v_congestion_cutoff = p_link->v_congestion_cutoff;
                        p_link->VDF_period[tau].v_congestion_cutoff = p_link->v_congestion_cutoff;

                        p_link->VDF_period[tau].Q_alpha = Q_alpha;
                        p_link->VDF_period[tau].Q_beta = Q_beta;
                        p_link->VDF_period[tau].Q_cd = Q_cd;
                        p_link->VDF_period[tau].Q_n = Q_n;
                        p_link->VDF_period[tau].Q_cp = Q_cp;
                        p_link->VDF_period[tau].Q_s = Q_s;
                    }
                    //we still update free speed for a comparison based on speed data
                    p_link->free_speed = highest_speed; // update or not, here is the question
                    p_link->VDF_period[tau].vf = highest_speed;
                    p_link->VDF_period[tau].v_congestion_cutoff = p_link->v_congestion_cutoff;

                    p_link->calculate_dynamic_VDFunction(0, false, p_link->vdf_type);


                    //            fprintf(p_file_tmc_link, "AM_t0,AM_t3,AM_P,AM_V,AM_D,AM_DC_ratio,AM_mu,AM_vu,AM_vf_reference, AM_v_mean,AM_vt2,AM_plf,");
                    float speed_reduction_factor = 0;
                    if (obs_P > 0.25 && t2_speed < p_link->v_congestion_cutoff)
                    {
                        speed_reduction_factor = p_link->v_congestion_cutoff / max(1.0f, t2_speed) - 1;
                    }
                    else
                    {
                        speed_reduction_factor = 0;
                    }

                    float queue_vc_delay_index = p_link->v_congestion_cutoff / max(1.0f, mean_speed_QVDF) - 1;
                    if (queue_vc_delay_index < 0)
                    {
                        queue_vc_delay_index = 0;
                    }

                    float BPR_vf_delay_index = max(highest_speed, (float)p_link->free_speed) / max(1.0f, mean_speed_BPR) - 1;
                    if (BPR_vf_delay_index < 0)
                    {
                        BPR_vf_delay_index = 0;
                    }

                    float log_VOC, log_DOC, log_P, log_sf, log_vfdi, log_vcdi;
                    log_VOC = log_DOC = log_P = log_sf = log_vfdi = log_vcdi = 0;

                    if (DOC_ratio > 0.00001)
                    {
                        log_DOC = log(DOC_ratio);
                    }

                    if (VOC_ratio > 0.00001)
                    {
                        log_VOC = log(VOC_ratio);
                    }

                    if (obs_P > 0.00001)
                    {
                        log_P = log(obs_P);
                    }


                    if (speed_reduction_factor > 0.00001)
                    {
                        log_sf = log(speed_reduction_factor);
                    }

                    if (BPR_vf_delay_index > 0.00001)
                    {
                        log_vfdi = log(BPR_vf_delay_index);
                    }

                    if (queue_vc_delay_index > 0.00001)
                    {
                        log_vcdi = log(queue_vc_delay_index);
                    }
                   // CBI
                    fprintf(p_file_tmc_link, "%f,%f,%f,%f,%f,",
                        obs_t0_in_hour, obs_t3_in_hour,  V, peak_hour_volume, D
                        );
                    fprintf(p_file_tmc_link, "%f,%f,%f,%f,%f,",
                        VOC_ratio, DOC_ratio, obs_P, speed_reduction_factor, BPR_vf_delay_index
                    );
                    fprintf(p_file_tmc_link, "%f,%f,%f,%f,%f,",
                        queue_vc_delay_index, mean_speed_BPR, mean_speed_QVDF, t2_speed, qdf

                    );
                    fprintf(p_file_tmc_link, "%f,%f,%f,%f,",
                        Q_n, Q_cp, Q_alpha, Q_beta
                    );

                    //QVDF and BPR model results
                    fprintf(p_file_tmc_link, "%f, %f,%f,%f,%f,%f, %f,%f,%f,%f,%f,",
                        g_link_vector[i].VDF_period[tau].link_volume,
                        
                        g_link_vector[i].VDF_period[tau].lane_based_D,
                        g_link_vector[i].VDF_period[tau].DOC,
                        g_link_vector[i].VDF_period[tau].P,
                        g_link_vector[i].VDF_period[tau].avg_queue_speed,
                        g_link_vector[i].VDF_period[tau].vt2,
                        
                        g_link_vector[i].VDF_period[tau].Q_mu,
                        g_link_vector[i].VDF_period[tau].Q_gamma,
                        g_link_vector[i].VDF_period[tau].lane_based_Vph,
                        g_link_vector[i].VDF_period[tau].VOC,
                        g_link_vector[i].VDF_period[tau].avg_speed_BPR
                    );

                }

                fprintf(p_file_tmc_link, "\"%s\",", g_link_vector[i].geometry.c_str());

                fprintf(p_file_tmc_link, "\"LINESTRING (");

                fprintf(p_file_tmc_link, "%f %f,", g_link_vector[i].TMC_from.x, g_link_vector[i].TMC_from.y);
                fprintf(p_file_tmc_link, "%f %f,", g_link_vector[i].TMC_to.x, g_link_vector[i].TMC_to.y);
                fprintf(p_file_tmc_link, ")\"");
                fprintf(p_file_tmc_link, ",");

                {
                    if (assignment.m_TMClink_map.find(g_link_vector[i].tmc_code) != assignment.m_TMClink_map.end())
                    {
                        int tmc_index = assignment.m_TMClink_map[g_link_vector[i].tmc_code];

                        double ObsSpeed[25];
                        int count_total = 0;
                        double EstSpeed[25], EstSpeedDiff[25];
                        double MAE_total = 0;
                        double MAPE_total = 0;
                        double RMSE_total = 0;


                        for (int t = 6 * 60; t < 20 * 60; t += 60)
                        {
                            int hour = t / 60;
                            ObsSpeed[hour] = g_TMC_vector[tmc_index].get_avg_hourly_speed(t);
                            fprintf(p_file_tmc_link, "%.1f,", ObsSpeed[hour]);
                        }

                        for (int t = 6 * 60; t < 20 * 60; t += 60)
                        {
                            double model_speed = g_link_vector[i].get_model_hourly_speed(t);
                            fprintf(p_file_tmc_link, "%.1f,", model_speed);
                            int hour = t / 60;
                            EstSpeed[hour] = model_speed;

                            if (EstSpeed[hour] > 1 && ObsSpeed[hour] >1 && analysis_hour_flag[hour] ==1)  // valid
                            {
                                EstSpeedDiff[hour] = fabs(EstSpeed[hour] - ObsSpeed[hour]);
                                MAE_total += fabs(EstSpeedDiff[hour]);
                                MAPE_total += fabs(EstSpeedDiff[hour]) / max(1.0, ObsSpeed[hour]);
                                RMSE_total += EstSpeedDiff[hour] * EstSpeedDiff[hour];
                                count_total += 1;
                            }
                            else
                                EstSpeedDiff[hour] = 0;

                        }

                        double MSE_total = RMSE_total / max(1, count_total);

                        fprintf(p_file_tmc_link, "%.2f,%.2f,%.2f,", MAE_total / max(1, count_total), MAPE_total / max(1, count_total) * 100, pow(MSE_total, 0.5));

                        for (int t = 6 * 60; t < 20 * 60; t += 60)
                        {
                            int hour = t / 60;
                            float volume = g_TMC_vector[tmc_index].get_avg_hourly_volume(t);
                            fprintf(p_file_tmc_link, "%.1f,", volume);
                        }
                        for (int t = 6 * 60; t < 20 * 60; t += 60)
                        {
                            int hour = t / 60;
                            ObsSpeed[hour] = g_TMC_vector[tmc_index].get_avg_hourly_speed(t);

                            float speed_ratio = ObsSpeed[hour] / max(1.0f, g_link_vector[i].TMC_highest_speed);

                            if (speed_ratio >= 1)
                                speed_ratio = 1;
                            fprintf(p_file_tmc_link, "%.3f,", speed_ratio);
                        }

                        for (int t = 6 * 60; t < 20 * 60; t += 5)
                        {

                            fprintf(p_file_tmc_link, "%.1f,", g_TMC_vector[tmc_index].get_avg_speed(t));
                        }

                        for (int t = 6 * 60; t < 20 * 60; t += 15)
                        {

                            fprintf(p_file_tmc_link, "%.1f,", g_TMC_vector[tmc_index].get_avg_speed(t));
                        }

                        for (int t = 6 * 60; t < 20 * 60; t += 5)
                        {
                            float speed = g_TMC_vector[tmc_index].get_avg_speed(t);
                            double volume = g_TMC_vector[tmc_index].get_avg_volume(t, &(g_link_vector[i]), speed, g_link_vector[i].TMC_highest_speed);
                            volume = volume * 12;  // compute hourly volume
                            fprintf(p_file_tmc_link, "%.2f,", volume);
                        }


                    }
                }
                fprintf(p_file_tmc_link, "\n");

            }

        }

        fclose(p_file_tmc_link);
    }
    else
    {
        dtalog.output() << "Error: File link_cbi_summary cannot be open." << endl;
        g_program_stop();
    }
}

void g_output_qvdf_file()
{
    // if 
    FILE* p_file_tmc_link = fopen("link_qvdf.csv", "w");


    if (p_file_tmc_link != NULL)
    {
        fprintf(p_file_tmc_link, "data_type,link_id,tmc_corridor_name,from_node_id,to_node_id,vdf_code,");

        for (int tau = 0; tau < min((size_t)3, assignment.g_DemandPeriodVector.size()); tau++)
        {
            fprintf(p_file_tmc_link, "QVDF_qdf%d,QVDF_n%d,QVDF_s%d,QVDF_cp%d,QVDF_cd%d,QVDF_alpha%d,QVDF_beta%d,", tau+1, tau + 1, tau + 1, tau + 1, tau + 1, tau + 1, tau + 1);
        }


        fprintf(p_file_tmc_link, "\n");

        std::map<__int64, int> TMC_long_id_mapping;  // this is used to mark if this cell_id has been identified or not

        //// sort data records
        for (int i = 0; i < g_link_vector.size(); i++)
        {
            if (g_link_vector[i].tmc_code.size() > 0)
            {
                __int64 TMC_long_key = (g_link_vector[i].tmc_corridor_id * 10000 + g_link_vector[i].tmc_road_sequence) * 10 + g_link_vector[i].link_seq_no;
                TMC_long_id_mapping[TMC_long_key] = g_link_vector[i].link_seq_no;
            }
        }


        std::map<__int64, int>::iterator it;

        for (it = TMC_long_id_mapping.begin(); it != TMC_long_id_mapping.end(); ++it)
        {
            int i = it->second;
            float highest_speed = g_link_vector[i].free_speed;

            if (assignment.m_TMClink_map.find(g_link_vector[i].tmc_code) != assignment.m_TMClink_map.end())
            {
                int tmc_index = assignment.m_TMClink_map[g_link_vector[i].tmc_code];

                if (g_TMC_vector[tmc_index].b_with_sensor_speed_data == false)
                    continue;

                highest_speed = g_TMC_vector[tmc_index].get_highest_speed();
            }
            else
            {
                continue; // no tmc record data
            }


            float free_speed = g_link_vector[i].free_speed;

            if (g_link_vector[i].lane_capacity > 5000)  // skip connectors
                continue;

            g_link_vector[i].update_kc(g_link_vector[i].free_speed);

//            fprintf(p_file_tmc_link, "link_id,tmc_corridor_name,from_node_id,to_node_id,vdf_code,");

            fprintf(p_file_tmc_link, "link,%s,%s,%d,%d,%s,",
                g_link_vector[i].link_id.c_str(),
                g_link_vector[i].tmc_corridor_name.c_str(),
                g_node_vector[g_link_vector[i].from_node_seq_no].node_id,
                g_node_vector[g_link_vector[i].to_node_seq_no].node_id,
                g_link_vector[i].vdf_code.c_str()
            );


           if (assignment.m_TMClink_map.find(g_link_vector[i].tmc_code) != assignment.m_TMClink_map.end())
            {
                int tmc_index = assignment.m_TMClink_map[g_link_vector[i].tmc_code];


                float highest_speed = 0;
                CLink* p_link = &(g_link_vector[i]);


                float updated_vc = g_TMC_vector[tmc_index].scan_highest_speed_and_vc(p_link->v_congestion_cutoff, p_link->free_speed, highest_speed);
                p_link->v_congestion_cutoff = updated_vc;
                p_link->update_kc(free_speed);

                for (int tau = 0; tau < min((size_t)3, assignment.g_DemandPeriodVector.size()); tau++)
                {

                    float assign_period_start_time_in_hour = assignment.g_DemandPeriodVector[tau].starting_time_slot_no * MIN_PER_TIMESLOT / 60.0;
                    float assign_period_end_time_in_hour = assignment.g_DemandPeriodVector[tau].ending_time_slot_no * MIN_PER_TIMESLOT / 60.0;
                    float assign_period_t2_peak_in_hour = assignment.g_DemandPeriodVector[tau].t2_peak_in_hour;

                    float obs_t0_in_hour, obs_t3_in_hour, obs_P;

                    obs_t0_in_hour = obs_t3_in_hour = obs_P = 0;

                    float V, D, VOC_ratio, DOC_ratio;
                    V = D = VOC_ratio = DOC_ratio = 0;

                    float peak_hour_volume = 0;
                    float mean_speed_BPR = 0;
                    float mean_speed_QVDF = 0;
                    float t2_speed = 0;
                    float plf = 1;
                    float qdf = 1;

                    float Q_n = 1;
                    float Q_s = 1;
                    float Q_cd = 1;
                    float Q_cp = 1;


                    obs_P = g_TMC_vector[tmc_index].scan_congestion_duration(tau,
                        assign_period_start_time_in_hour,
                        assign_period_end_time_in_hour,
                        assign_period_t2_peak_in_hour,
                        p_link->v_congestion_cutoff, p_link,
                        obs_t0_in_hour, obs_t3_in_hour, obs_P,
                        V, peak_hour_volume, D, VOC_ratio, DOC_ratio,
                        mean_speed_BPR, mean_speed_QVDF, highest_speed, t2_speed,
                        plf, qdf,
                        Q_n, Q_s, Q_cd, Q_cp);

                    qdf = plf;

                    double Q_alpha = 8.0 / 15 * Q_cp * pow(Q_cd, Q_s);
                    double Q_beta = Q_n * Q_s;

                    Q_alpha = g_TMC_vector[tmc_index].check_feasible_range(Q_alpha, 0.27, 0.01, 1);
                    Q_beta = g_TMC_vector[tmc_index].check_feasible_range(Q_beta, 1.14, 0.5, 5);

                    p_link->VDF_period[tau].Q_alpha = Q_alpha;
                    p_link->VDF_period[tau].Q_beta = Q_beta;
                    p_link->VDF_period[tau].Q_cd = Q_cd;
                    p_link->VDF_period[tau].Q_n = Q_n;
                    p_link->VDF_period[tau].Q_cp = Q_cp;
                    p_link->VDF_period[tau].Q_s = Q_s;

                    if (p_link->vdf_code.size() > 0)
                    {
                        g_vdf_type_map[p_link->vdf_code].record_qvdf_data(p_link->VDF_period[tau], tau);
                    }

                        g_vdf_type_map["all"].record_qvdf_data(p_link->VDF_period[tau], tau);

//                    fprintf(p_file_tmc_link, "QVDF_qdf%d,QVDF_n%d,QVDF_s%d,QVDF_cp%d,QVDF_cd%d,QVDF_alpha%d,QVDF_beta%d\n");
                    fprintf(p_file_tmc_link, "%f,%f,%f,%f,%f,%f,%f,",
                        qdf, Q_n, Q_s, Q_cp, Q_cd, Q_alpha, Q_beta 
                    );

                }


                fprintf(p_file_tmc_link, "\n");

            }

        }

        {
            std::map<string, CVDF_Type>::iterator it;
            for (it = g_vdf_type_map.begin(); it != g_vdf_type_map.end(); ++it)
            {
                fprintf(p_file_tmc_link, "vdf_code,,,,,%s,", it->first.c_str());
                for (int tau = 0; tau < min((size_t)3, assignment.g_DemandPeriodVector.size()); tau++)
                {
                    it->second.computer_avg_parameter(tau);
                    fprintf(p_file_tmc_link, "%f,%f,%f,%f,%f,%f,%f,",
                        it->second.VDF_period_sum[tau].queue_demand_factor,
                        it->second.VDF_period_sum[tau].Q_n,
                        it->second.VDF_period_sum[tau].Q_s,
                        it->second.VDF_period_sum[tau].Q_cp,
                        it->second.VDF_period_sum[tau].Q_cd,
                        it->second.VDF_period_sum[tau].Q_alpha,
                        it->second.VDF_period_sum[tau].Q_beta);

                }
                fprintf(p_file_tmc_link, "\n");

            }
        }

        fclose(p_file_tmc_link);
    }
    else
    {
        dtalog.output() << "Error: File link_qvdf cannot be open." << endl;
        g_program_stop();

    }

}