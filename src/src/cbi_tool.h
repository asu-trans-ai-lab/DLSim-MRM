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

class CTMC_Link
{
public:

    CTMC_Link()
    {
        b_with_sensor_speed_data = false;
        for (int t = 0; t < MAX_TIMEINTERVAL_PerDay; ++t)
        {
            speed_sum[t] = -1;
            avg_speed[t] = -1;
            est_speed[t] = -1;
            speed_lowest[t] = 99999;
            count[t] = 0;
        }


        //for (int d = 0; d < MAX_DAY_PerYear; d++)
        //{
        //    for (int t = 0; t < MAX_TIMEINTERVAL_PerDay; t++)
        //    {
        //        speed_day[d][t] = -1;
        //    }
        //}

    }

    ~CTMC_Link()
    {
    }

    float speed_lowest[MAX_TIMEINTERVAL_PerDay];
    float volume_per_hour_per_lane[MAX_TIMEINTERVAL_PerDay];
    float count[MAX_TIMEINTERVAL_PerDay];


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


    float scan_congestion_duration(int peak_no, int starting_time_in_hour, int ending_time_in_hour, float& FD_vc, float& obs_t0_in_hour, float& obs_t3_in_hour,
        CLink* pLink, float& assignment_volume, float& assignment_VMT, float& assignment_VHT, float& assignment_VDT, float& assignment_VCDT, float& congestion_demand, float& congestion_dc_ratio, float& mean_congestion_flow_rate_mu, float& mean_congestion_speed,
        float& highest_speed, float& mean_speed, float& t2_speed, float& t2_queue, float& gamma)
    {

        obs_t0_in_hour = -1;
        obs_t3_in_hour = -1;
        est_PLF[peak_no] = 1;
        est_QDF_n[peak_no] = 1.0804;
        est_QDF_s[peak_no] = 0.4743;
        est_QDF_beta[peak_no] = 1.0804*0.4743;
        est_gamma[peak_no] = 0;

        assignment_VHT = 0;
        assignment_VMT = 0;
        assignment_VDT = 0;

        int obs_t0_in_interval = starting_time_in_hour * 12;
        int obs_t3_in_interval = ending_time_in_hour * 12;

        float total_speed_value = 0;
        float total_speed_count = 0;

        highest_speed = 0;


        for (int t_in_min = 6 * 60; t_in_min < 20 * 60; t_in_min += 5)
        {

            float avg_speed = record_avg_speed(t_in_min);

            if (avg_speed > highest_speed)
                highest_speed = avg_speed;
        }


        for (int t_in_min = starting_time_in_hour * 60; t_in_min < ending_time_in_hour * 60; t_in_min += 5)
        {

            float avg_speed = record_avg_speed(t_in_min);

            total_speed_value += avg_speed;
            total_speed_count++;
        }


        // use the measured speed to update 
        pLink->update_kc(highest_speed);
        FD_vc = pLink->vc;

        mean_speed = total_speed_value / max(1, total_speed_count);

        int t_mid = (starting_time_in_hour + ending_time_in_hour) * 12 / 2;
        est_PLF[peak_no] = max(0.1, ending_time_in_hour - starting_time_in_hour); // number of hours in peak period
        assignment_volume = 0;

        int t_lowest_speed = -1;
        double lowest_speed = 99999;


        if (avg_speed[t_mid] > 0)
        {
            float VHT_sum = 0;
            float VMT_sum = 0;
            float VDT_sum = 0;
            float VCDT_sum = 0; // congestion_delay
            for (int t_in_min = starting_time_in_hour * 60; t_in_min < ending_time_in_hour * 60; t_in_min += 5)
            {

                float avg_speed = record_avg_speed(t_in_min);
                float volume = pLink->get_volume_from_speed(avg_speed, highest_speed);
                float VHT = volume / 12 * pLink->link_distance_in_km / avg_speed;
                float VMT = volume / 12 * pLink->link_distance_in_km;
                float VDT = volume / 12 * max(0, pLink->link_distance_in_km / avg_speed - pLink->link_distance_in_km / (max(1, highest_speed)));  //base is free-flow travel time 
                float VCDT = volume / 12 * max(0, pLink->link_distance_in_km / avg_speed - pLink->link_distance_in_km / (max(1, pLink->vc)));  // base is vc; speed at capacity
                VHT_sum += VHT * pLink->number_of_lanes;
                VMT_sum += VMT * pLink->number_of_lanes;
                VDT_sum += VDT * pLink->number_of_lanes;
                VCDT_sum += VCDT * pLink->number_of_lanes;
                assignment_volume += volume / 12; // 12 5-min interval per hour


                if (avg_speed < lowest_speed)
                {
                    t_lowest_speed = t_in_min / 5;
                    lowest_speed = avg_speed;
                }



            }

            assignment_VHT = VHT_sum;
            assignment_VMT = VMT_sum;
            assignment_VDT = VDT_sum;
            assignment_VCDT = VCDT_sum;
            est_vcd[peak_no] = VCDT_sum;

        }

        if (avg_speed[t_mid] < 0)
            return 0;

        if (avg_speed[t_mid] >= FD_vc)
        {
            // condition 1: the congestion in the middle 
            //not met
            // condition 2: 
            if (lowest_speed <= FD_vc)
            {
                t_mid = t_lowest_speed; // update t_mid using lowest speed time index.
            }
            else
            {
                return 0;
            }
        }


        // below is congested period analysis

        bool b_t0_found_flag = false;
        bool b_t3_found_flag = false;

        if (avg_speed[t_mid] < FD_vc)  // exit if the speed is higher than vc
        {
            obs_t3_in_interval = t_mid;
            obs_t3_in_hour = t_mid * 1.0 / 12.0;
            b_t3_found_flag = true;

            obs_t0_in_interval = t_mid;
            b_t0_found_flag = true;
            obs_t0_in_hour = t_mid * 1.0 / 12.0;

        }

        for (int t = t_mid + 1; t < 24 * 12; t += 1)  // move forward from mid t
        {
            if (avg_speed[t] < 1)
            {
                int i_no_data = 1;

                //                dtalog.output() << " Error: TMC_link " << pLink->tmc_code.c_str() << " has no data at timt t " << t << g_time_coding(t * 5).c_str() << endl;;
                break;
            }

            obs_t3_in_interval = t;
            obs_t3_in_hour = t * 1.0 / 12.0;
            b_t3_found_flag = true;

            if (avg_speed[t] > FD_vc)  // exit if the speed is higher than vc
            {
                break;
            }
        }


        for (int t = t_mid - 1; t >= 0 * 12; t -= 1) // move backward from t_mid 
        {
            if (avg_speed[t] < 1)
            {
                int i_no_data = 1;
                dtalog.output() << " Error: TMC link " << pLink->tmc_code.c_str() << " has no data at timt t =" << t * 5 << " min" << endl;;
                break;
            }
            if (avg_speed[t] > FD_vc)
            {

                break;
            }

            obs_t0_in_interval = t;
            b_t0_found_flag = true;
            obs_t0_in_hour = t * 1.0 / 12.0;

        }


        if (b_t0_found_flag == false || b_t3_found_flag == false)  // two boundaries are found
            return 0;


        congestion_demand = 0;
        total_speed_count = 0;  // initial values
        total_speed_value = 0;

        for (int t = obs_t0_in_interval; t <= obs_t3_in_interval; t += 1) // move between congestion duration per interval
        {
            if (avg_speed[t] >= 1)
            {
                total_speed_value += avg_speed[t];
                float volume = pLink->get_volume_from_speed(avg_speed[t], highest_speed);
                congestion_demand += volume / 12; // 12 5-min interval per hour
                total_speed_count++;
            }

        }

        // test
        float test_volume = pLink->get_volume_from_speed(avg_speed[t_mid], highest_speed);

        mean_congestion_flow_rate_mu = congestion_demand / max(1, total_speed_count) * 12;
        float    obs_P_in_hour = (obs_t3_in_hour - obs_t0_in_hour);  // congestion duration P

        congestion_demand = mean_congestion_flow_rate_mu * obs_P_in_hour;  // recalculate D in case there are missing data. 

        congestion_dc_ratio = congestion_demand / max(1, pLink->lane_capacity);  // unit: demand: # of vehicles, lane_capacity # of vehicles per hours: dc ratio has a unit of hour, but it is different from P

        mean_congestion_speed = total_speed_value / max(1, total_speed_count);

        t2_speed = avg_speed[t_mid];  // if we use a pure second order model, we should consider t2= 2/3(t3-t0)+ t0
        float t2_waiting_time = 1 / t2_speed - 1 / FD_vc;  // assume 1 mile road, unit: 1 hour
        t2_queue = t2_waiting_time * mean_congestion_flow_rate_mu;
        gamma = t2_queue * 12 / pow((obs_t3_in_hour - obs_t0_in_hour) / 2.0, 4);


        // calibration
        if (obs_P_in_hour > 5) // setup upper bound
            obs_P_in_hour = 5;

        if (obs_P_in_hour > 1)  // only for congested hour
        {
            est_PLF[peak_no] = max(1, assignment_volume / max(pLink->lane_capacity, congestion_demand));
            est_QDF_n[peak_no] = max(1, log(obs_P_in_hour) / log(congestion_dc_ratio));
            float Vmin_ratio = FD_vc / max(1, t2_speed);
            est_QDF_s[peak_no] = max(0.0001, log(Vmin_ratio) / log(obs_P_in_hour));
            est_QDF_beta[peak_no] = est_QDF_s[peak_no] * est_QDF_s[peak_no];
            est_gamma[peak_no] = gamma;
        }


        return obs_P_in_hour;
    }

    string tmc_code;

    int matching_link_no;
    bool b_with_sensor_speed_data;

    // construction
    void add_speed_data(int day_no, int time_in_min, float speed_value)
    {
        b_with_sensor_speed_data = true;
        int t = time_in_min / 5;

        if (t >= 0 && t < MAX_TIMEINTERVAL_PerDay)
        {
            //            speed_day[day_no][t] = speed_value;

            if (speed_sum[t] < 0 || count[t] == 0)
                speed_sum[t] = 0;

            speed_sum[t] += speed_value;
            count[t] += 1;

            if (speed_value < speed_lowest[t])
                speed_lowest[t] = speed_value;

        }

    }

    float record_avg_speed(int time_in_min)
    {
        int t = time_in_min / 5;

        if (count[t] == 0)
            return -1;

        if (t >= 0 && t < MAX_TIMEINTERVAL_PerDay)
        {

            avg_speed[t] = speed_sum[t] / max(1, count[t]);
            return avg_speed[t];
        }

        return -1;
    }

    float get_avg_speed(int time_in_min)
    {
        int t = time_in_min / 5;
        if (t >= 0 && t < MAX_TIMEINTERVAL_PerDay)
        {

            return speed_sum[t] / max(1, count[t]);

        }
        return -1;
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
                total_speed_value += speed_sum[t + tt] / max(1, count[t + tt]);
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
                total_speed_value += speed_sum[t + tt] / max(1, count[t + tt]);
                total_speed_count++;
            }
        }

        return total_speed_value / max(1, total_speed_count);
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
    //float speed_day[MAX_DAY_PerYear][MAX_TIMEINTERVAL_PerDay];

    float est_speed[MAX_TIMEINTERVAL_PerDay];
    float est_volume_per_hour_per_lane[MAX_TIMEINTERVAL_PerDay];

    float pred_speed[MAX_TIMEINTERVAL_PerDay];
    float pred_per_hour_per_lane[MAX_TIMEINTERVAL_PerDay];

    // AM 0: noon 1, PM: 2

    float est_PLF[MAX_TIMEPERIODS];
    float est_QDF_n[MAX_TIMEPERIODS];  // P= (D/C)^n
    float est_QDF_s[MAX_TIMEPERIODS];  // vc/vt2= (P)^s, P Vmin slope
    float est_QDF_beta[MAX_TIMEPERIODS];  
    // first layer of varaibles 
    float est_V[MAX_TIMEPERIODS];
    float est_D[MAX_TIMEPERIODS];
    float est_DCRatio[MAX_TIMEPERIODS];

    // 2nd layer of varaibles 
    float est_P[MAX_TIMEPERIODS];
    float est_mu[MAX_TIMEPERIODS];
    float est_vcd[MAX_TIMEPERIODS];
    float est_vt2[MAX_TIMEPERIODS];

    float est_gamma[MAX_TIMEPERIODS];
    float est_t0[MAX_TIMEPERIODS];
    float est_t3[MAX_TIMEPERIODS];

    float get_est_speed(int time_in_min)
    {
        int t = time_in_min / 5;
        if (t >= 0 && t < MAX_TIMEINTERVAL_PerDay)
        {
            if (est_speed[t] >= 1)
                return est_speed[t];
        }
        return 0;
    }
    float get_est_hourly_speed(int time_in_min)
    {
        int t = time_in_min / 5;
        float total_speed_value = 0;
        int total_speed_count = 0;

        for (int tt = 0; tt < 12; tt++)
        {

            if (t + tt >= 0 && t + tt < MAX_TIMEINTERVAL_PerDay)
            {
                if (est_speed[t + tt] >= 1)
                {
                    total_speed_value += est_speed[t + tt];
                    total_speed_count++;
                }
            }
        }

        return total_speed_value / max(1, total_speed_count);
    }
    float get_pred_speed(int time_in_min)
    {
        int t = time_in_min / 5;
        if (t >= 0 && t < MAX_TIMEINTERVAL_PerDay)
        {
            if (pred_speed[t] >= 1)
                return pred_speed[t];
        }
        return 0;
    }
    float get_pred_hourly_speed(int time_in_min)
    {
        int t = time_in_min / 5;
        float total_speed_value = 0;
        int total_speed_count = 0;

        for (int tt = 0; tt < 12; tt++)
        {

            if (t + tt >= 0 && t + tt < MAX_TIMEINTERVAL_PerDay)
            {
                if (pred_speed[t + tt] >= 1)
                {
                    total_speed_value += pred_speed[t + tt];
                    total_speed_count++;
                }
            }
        }

        return total_speed_value / max(1, total_speed_count);
    }

    float perform_estimation(bool bPredictionFlag, CLink* pLink, int peak_no, float assign_period_start_time_in_hour, float assign_period_end_time_in_hour, float t2, float V, float laneCapacity, float vf, float vcd, float vct, float& DTASpeed, float& DTAP, float& TDAvgSpeedDiff)
    {
        int p = peak_no;

        TDAvgSpeedDiff = -1;
        est_V[p] = V;
        est_D[p] = V / est_PLF[p];
        double coef_a = pow(pLink->kc, pLink->s3_m);
        double coef_b = pow(pLink->kc, pLink->s3_m) * pow(vf, pLink->s3_m / 2.0);
        double coef_c = pow(est_D[p], pLink->s3_m);

        DTASpeed = pow((coef_b + pow(coef_b * coef_b - 4.0 * coef_a * coef_c, 0.5)) / (2.0 * coef_a), 2.0 / pLink->s3_m);    //under uncongested condition

        est_DCRatio[p] = est_D[p] / laneCapacity;
        est_P[p] = pow(est_DCRatio[p], est_QDF_n[p]);
        est_mu[p] = pow(est_DCRatio[p], est_QDF_n[p] * (-1)) * est_D[p];
        est_vt2[p] = vcd / pow(est_P[p], est_QDF_s[p]); //est_vt2=vc/P^s   

        DTAP = est_P[p];
        // setup up default values
        for (int t_in_min = assign_period_start_time_in_hour * 60; t_in_min < assign_period_end_time_in_hour * 60; t_in_min += 5)
        {
            int t_interval = t_in_min / 5;

            if (bPredictionFlag)
                pred_speed[t_interval] = DTASpeed;
            else
                est_speed[t_interval] = DTASpeed;
        }

        if (DTAP < 1) // for P < 1, we will skip the creation of curve below vc. 
            return DTASpeed;

        float est_q_t2 = 1.0 / est_vt2[p] * est_mu[p];
        est_gamma[p] = est_q_t2 * 4 / pow(est_P[p] / 2, 4);  // because q_tw = w*mu =1/4 * gamma (P/2)^4, => 1/vt2 * mu = 1/4 * gamma  * (P/2)^4

        est_t0[p] = t2 - 0.5 * est_P[p];
        est_t3[p] = t2 + 0.5 * est_P[p];;

        //if (est_P[p] > 2 && est_DCRatio[p] > 2)
        //{
        //    int idebug = 1;
        //}


        for (int t_in_min = assign_period_start_time_in_hour * 60; t_in_min < assign_period_end_time_in_hour * 60; t_in_min += 5)
        {
            float t = t_in_min / 60.0;  // t in hour
            float td_queue;

            if (est_t0[p] <= t && t <= est_t3[p])
                td_queue = 0.25 * est_gamma[p] * pow((t - est_t0[p]), 2) * pow((t - est_t3[p]), 2);
            else
                td_queue = 0;

            float td_speed;
            if (est_P[p] < 0.25)
            {
                if (t < est_t0[p])
                {
                    td_speed = vf - ((vf - vct) / (est_t0[p] - assign_period_start_time_in_hour)) * (t - assign_period_start_time_in_hour);
                }
                else if (t > est_t3[p])
                {
                    td_speed = vf - ((vf - vct) / (assign_period_end_time_in_hour - est_t3[p])) * (assign_period_end_time_in_hour - t);
                }
                else
                    td_speed = vcd;
            }
            else // est_P > 0.25
            {
                if (t < est_t0[p])
                    td_speed = vf - ((vf - vcd) / (est_t0[p] - assign_period_start_time_in_hour)) * (t - assign_period_start_time_in_hour);
                else if (t > est_t3[p])
                    td_speed = vf - ((vf - vcd) / (assign_period_end_time_in_hour - est_t3[p])) * (assign_period_end_time_in_hour - t);
                else

                    td_speed = 1 / ((td_queue / est_mu[p]) + (1 / vcd));

                //                cout << "td_queue t" << t << " =  " << td_queue << ", speed =" << td_speed << endl;
            }

            int t_interval = t_in_min / 5;

            if (bPredictionFlag)
                pred_speed[t_interval] = td_speed;
            else
                est_speed[t_interval] = td_speed;

        }

        if (bPredictionFlag == false)
        {
            // tally the mean DTA speed
            float est_speed_total = 0;
            int est_speed_count = 0;
            float total_speed_abs_diff = 0;
            for (int t_in_min = assign_period_start_time_in_hour * 60; t_in_min < assign_period_end_time_in_hour * 60; t_in_min += 5)
            {
                int t_interval = t_in_min / 5;
                est_speed_total += est_speed[t_interval];
                est_speed_count++;

                if (avg_speed[t_interval] >= 1)  // with data
                {
                    total_speed_abs_diff += fabs(est_speed[t_interval] - avg_speed[t_interval]);
                }
            }

            TDAvgSpeedDiff = total_speed_abs_diff / max(1, est_speed_count);
            DTASpeed = est_speed_total / max(1, est_speed_count);
            return est_vt2[p];
        }
        else
        {
            // tally the mean DTA speed
            float pred_speed_total = 0;
            int pred_speed_count = 0;
            float total_speed_abs_diff = 0;
            for (int t_in_min = assign_period_start_time_in_hour * 60; t_in_min < assign_period_end_time_in_hour * 60; t_in_min += 5)
            {
                int t_interval = t_in_min / 5;
                pred_speed_total += pred_speed[t_interval];
                pred_speed_count++;

                if (est_speed[t_interval] >= 1)  // with data
                {
                    total_speed_abs_diff += fabs(pred_speed[t_interval] - est_speed[t_interval]);
                }
            }

            TDAvgSpeedDiff = total_speed_abs_diff / max(1, pred_speed_count);
            DTASpeed = pred_speed_total / max(1, pred_speed_count);
            return est_vt2[p];

        }

    }


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

    int char_link_distance_in_km = strlen(string_line);

    char ch, buf_hh[32] = { 0 }, buf_mm[32] = { 0 }, buf_ss[32] = { 0 };

    int yyyy, month, day;
    char dd1, dd2, hh1, hh2, mm1, mm2;
    float ddf1, ddf2, hhf1, hhf2, mmf1, mmf2;
    float global_minute = 0;
    float dd = 0, hh = 0, mm = 0;
    int i = 11;  //2019-01-01 06:00:00,
    int buffer_i = 0, buffer_k = 0, buffer_j = 0;
    int num_of_colons = 0;
    int num_of_underscore = 0;

    int char_link_distance_in_km_yymmdd = 10;
    i = 0;

    //    sscanf(string_line, "%d/%d/%d", &month, &day, &yyyy );
    sscanf(string_line, "%d-%d-%d", &yyyy, &month, &day);



    day_of_week_flag = g_dayofweek(yyyy, month, day);
    day_of_year = g_dayofyear(yyyy, month, day);

    g_dayDataMap[day_of_year] = month * 100 + day;


    /// 
    i = 11;

    while (i < char_link_distance_in_km)
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

        if (i == char_link_distance_in_km) //start a new time string
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
            float speed = 0;
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

            parser_measurement.GetValueByFieldName("speed", speed, false);

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

            g_TMC_vector[index].add_speed_data(day_of_year, global_time, speed);



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

    }

    return false;

}
void g_output_tmc_file(bool bReadingDataReady)
{
    FILE* p_file_tmc_link = fopen("link_cbi_summary.csv", "w");
    FILE* p_file_tmc_vdf = fopen("link_cbi_vdf.csv", "w");

    if (p_file_tmc_vdf == NULL)
    {
        dtalog.output() << "Error: File link_cbi_vdf.csv cannot be opened.\n It might be currently used and locked by EXCEL." << endl;
        g_program_stop();
        return;
    }

    if (p_file_tmc_link != NULL)
    {
        fprintf(p_file_tmc_link, "link_id,tmc,tmc_corridor_name,tmc_corridor_id,tmc_road_order,tmc_road_sequence,tmc_road,tmc_direction,tmc_intersection,tmc_highest_speed,link_no,from_node_id,to_node_id,link_type,");
        fprintf(p_file_tmc_link, "link_type_code, FT, AT, nlanes,link_distance_in_km,free_speed,capacity,kc, ");
        fprintf(p_file_tmc_link, "AM_vc,AM_QHF,AM_QDF_n,AM_QDF_s,AM_t0, AM_t3, AM_P, AM_Assign_V,AM_Assign_VMT,AM_Assign_VHT, AM_Assign_VDT, AM_Assign_VCDT,AM_D, AM_DC_ratio, AM_mu, AM_vu, AM_vf_reference, AM_v_mean, AM_t2_v, AM_t2_queue, AM_gamma, AM_DTASpeed1, AM_DTAP1, AM_DTATDSpdDiff,");
        fprintf(p_file_tmc_link, "MD_vc,MD_QHF,MD_QDF_n,MD_QDF_s,MD_t0, MD_t3, MD_P, MD_Assign_V,MD_Assign_VMT,MD_Assign_VHT, MD_Assign_VDT, MD_Assign_VCDT,MD_D, MD_DC_ratio, MD_mu, MD_vu, MD_vf_reference, MD_v_mean, MD_t2_v, MD_t2_queue, MD_gamma, MD_DTASpeed1, MD_DTAP1, MD_DTATDSpdDiff,");
        fprintf(p_file_tmc_link, "PM_vc,PM_QHF,PM_QDF_n,PM_QDF_s,PM_t0, PM_t3, PM_P, PM_Assign_V,PM_Assign_VMT,PM_Assign_VHT, PM_Assign_VDT, PM_Assign_VCDT,PM_D, PM_DC_ratio, PM_mu, PM_vu, PM_vf_reference, PM_v_mean, PM_t2_v, PM_t2_queue, PM_gamma,");
        fprintf(p_file_tmc_link, "geometry, tmc_geometry, STA_speed1,STA_VOC1,STA_volume1,STA_speed2,STA_VOC2,STA_volume2,STA_speed3,STA_VOC3,STA_volume3,STA_speed4,STA_VOC4,STA_volume4,");

        fprintf(p_file_tmc_vdf, "link_id,from_node_id,to_node_id,link_type,free_speed,capacity,kc,");

        for (int k = 1; k <= 3; k++)
        {
            fprintf(p_file_tmc_vdf, "VDF_plf%d,VDF_m%d,VDF_n%d,VDF_s%d,VDF_beta%d,VDF_gamma%d,VDF_sav%d,QVDF_d%d,QVDF_dc%d,QVDF_p%d,QVDF_mu%d,QVDF_vcd%d,",
                                    k,  k, k,        k,        k,          k,     k,k,    k,     k,       k,      k );
        }

        fprintf(p_file_tmc_vdf, "\n");

   
        //if(bReadingDataReady)
        {
            for (int t = 6 * 60; t < 20 * 60; t += 60)
            {
                int hour = t / 60;
                int minute = t - hour * 60;
                fprintf(p_file_tmc_link, "vh%02d,", hour, minute);
            }
            for (int t = 6 * 60; t < 20 * 60; t += 60)
            {
                int hour = t / 60;
                int minute = t - hour * 60;

                fprintf(p_file_tmc_link, "evh%02d,", hour);
            }

            for (int t = 6 * 60; t < 20 * 60; t += 60)
            {
                int hour = t / 60;
                int minute = t - hour * 60;

                fprintf(p_file_tmc_link, "evh%02ddiff,", hour);
            }

            fprintf(p_file_tmc_link, "evhMAE,evhMAPE,evhRMSE,");

            for (int t = 6 * 60; t < 20 * 60; t += 60)
            {
                int hour = t / 60;
                int minute = t - hour * 60;
                fprintf(p_file_tmc_link, "vhr%02d,", hour, minute);
            }
            for (int t = 6 * 60; t < 20 * 60; t += 60)
            {
                int hour = t / 60;
                int minute = t - hour * 60;

                fprintf(p_file_tmc_link, "evhr%02d,", hour);
            }
            fprintf(p_file_tmc_link, "STA_AM_AE,STA_MD_AE,STA_PM_AE,STA_MAE,STA_MAPE,STA_RMSE,");


            for (int t = 6 * 60; t < 20 * 60; t += 5)
            {
                int hour = t / 60;
                int minute = t - hour * 60;

                fprintf(p_file_tmc_link, "v%02d:%02d,", hour, minute);
            }

            for (int t = 0 * 60; t < 24 * 60; t += 15)
            {
                int hour = t / 60;
                int minute = t - hour * 60;

                fprintf(p_file_tmc_link, "v%02d:%02d,", hour, minute);
            }

            for (int t = 6 * 60; t < 20 * 60; t += 15)
            {
                int hour = t / 60;
                int minute = t - hour * 60;

                fprintf(p_file_tmc_link, "vr%02d:%02d,", hour, minute);
            }


            for (int t = 6 * 60; t < 20 * 60; t += 15)
            {
                int hour = t / 60;
                int minute = t - hour * 60;

                fprintf(p_file_tmc_link, "q%02d:%02d,", hour, minute);
            }
            for (int t = 6 * 60; t < 20 * 60; t += 15)
            {
                int hour = t / 60;
                int minute = t - hour * 60;

                fprintf(p_file_tmc_link, "q%02d:%02d,", hour, minute);
            }
            /// <summary>
            /// / etimation 
            /// </summary>
            for (int t = 6 * 60; t < 20 * 60; t += 15)
            {
                int hour = t / 60;
                int minute = t - hour * 60;

                fprintf(p_file_tmc_link, "ev%02d:%02d,", hour, minute);
            }
        }


        fprintf(p_file_tmc_link, "\n");

        std::map<_int64, int> TMC_long_id_mapping;  // this is used to mark if this cell_id has been identified or not

        //// sort data records
        for (int i = 0; i < g_link_vector.size(); i++)
        {
            if (g_link_vector[i].link_id == "5666")
            {
                int idebug = 1;
            }

            if (g_link_vector[i].tmc_code.size() > 0)
            {

                _int64 TMC_long_key = (g_link_vector[i].tmc_corridor_id * 10000 + g_link_vector[i].tmc_road_sequence) * 10 + g_link_vector[i].link_seq_no;
                TMC_long_id_mapping[TMC_long_key] = g_link_vector[i].link_seq_no;
            }
        }


        std::map<_int64, int>::iterator it;

        for (it = TMC_long_id_mapping.begin(); it != TMC_long_id_mapping.end(); ++it)
        {
            int i = it->second;
            float highest_speed = g_link_vector[i].free_speed;

            if (assignment.m_TMClink_map.find(g_link_vector[i].tmc_code) != assignment.m_TMClink_map.end())
            {
                int tmc_index = assignment.m_TMClink_map[g_link_vector[i].tmc_code];

                highest_speed = g_TMC_vector[tmc_index].get_highest_speed();
                //    if (g_TMC_vector[tmc_index].b_with_sensor_speed_data == false)  // with reading speed data only
                //        continue;
            }


            float free_speed = g_link_vector[i].free_speed;
            g_link_vector[i].update_kc(g_link_vector[i].free_speed);


            fprintf(p_file_tmc_link, "%s,%s,%s,%d,%d,%d,%s,%s,%s,%.2f,%d,%d,%d,%d,%s,%d,%d,%d,%f,%f,%f,%f,",
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
                g_link_vector[i].number_of_lanes,
                g_link_vector[i].link_distance_in_km,
                g_link_vector[i].free_speed,
                g_link_vector[i].lane_capacity,
                g_link_vector[i].kc

            );


            if (assignment.m_TMClink_map.find(g_link_vector[i].tmc_code) == assignment.m_TMClink_map.end())
            {
                fprintf(p_file_tmc_link, "\n");
                continue;

            }
            float FD_vc = 0;
            float obs_t0_in_hour = -1;
            float obs_t3_in_hour = -1;
            float mean_congestion_speed = 0;
            float congestion_dc_ratio = 0;
            float obs_P_in_hour = 0;
            float mean_congestion_mu = 0;
            highest_speed = 0;
            float mean_speed = 0;
            float congestion_D = 0;
            float assignment_D = 0;
            float assignment_VHT = 0;
            float assignment_VMT = 0;
            float assignment_VDT = 0;
            float assignment_VCDT = 0;

            float t2_speed = 0;
            float t2_queue = 0;
            float gamma = 0;
            float assign_period_start_time_in_hour = 6;
            float assign_period_end_time_in_hour = 11;
            int starting_time_in_hour = 6;
            int ending_time_in_hour = 10;
            float t2 = 7.5;
            float DTASpeed = -1;
            float DTAP = -1;
            float TDAvgSpeedDiff = -1;

            float STA_AM_MAE = -1;
            float STA_MD_MAE = -1;
            float STA_PM_MAE = -1;
            float STA_AM_APE = -1;
            float STA_MD_APE = -1;
            float STA_PM_APE = -1;
            // P D analysis


            if (assignment.m_TMClink_map.find(g_link_vector[i].tmc_code) != assignment.m_TMClink_map.end())
            {
                int tmc_index = assignment.m_TMClink_map[g_link_vector[i].tmc_code];

                if (g_link_vector[i].tmc_code == "115-04402")
                {
                    int idebug = 1;
                }
                int peak_no = 0;  //AM

                assign_period_start_time_in_hour = 6;
                assign_period_end_time_in_hour = 10;
                starting_time_in_hour = 6;
                ending_time_in_hour = 9;


                CLink* pLink = &(g_link_vector[i]);

                if (pLink->link_id == "213")
                    int debug_flag = 1;

                obs_P_in_hour = g_TMC_vector[tmc_index].scan_congestion_duration(peak_no, starting_time_in_hour,
                    ending_time_in_hour, FD_vc, obs_t0_in_hour, obs_t3_in_hour, pLink, assignment_D, assignment_VMT, assignment_VHT, assignment_VDT, assignment_VCDT, congestion_D,
                    congestion_dc_ratio, mean_congestion_mu, mean_congestion_speed, highest_speed, mean_speed, t2_speed, t2_queue, gamma);

                float V = assignment_D;
                float laneCapacity = pLink->lane_capacity;
                float vf = highest_speed;
                float vcd = pLink->vc;
                float vct = pLink->vc;
                t2 = 7.75;  // calibrated based on 101 data

                g_TMC_vector[tmc_index].perform_estimation(false, pLink, peak_no, assign_period_start_time_in_hour, assign_period_end_time_in_hour, t2, V, laneCapacity, vf, vcd, vct, DTASpeed, DTAP, TDAvgSpeedDiff);

                STA_AM_MAE = fabs(g_link_vector[i].VDF_STA_speed[1] - mean_speed);
                STA_AM_APE = STA_AM_MAE / max(1, mean_speed);

                float QHF = 3;
                if (congestion_D > 1 && obs_P_in_hour >= 1)
                    QHF = V / congestion_D;

                fprintf(p_file_tmc_link, "%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,",
                    FD_vc,
                    QHF, g_TMC_vector[tmc_index].est_QDF_n[peak_no], g_TMC_vector[tmc_index].est_QDF_s[peak_no],
                    obs_t0_in_hour, obs_t3_in_hour, obs_P_in_hour, assignment_D, assignment_VMT, assignment_VHT, assignment_VDT, assignment_VCDT,
                    congestion_D, congestion_dc_ratio, mean_congestion_mu, mean_congestion_speed,
                    highest_speed, mean_speed,
                    t2_speed, t2_queue, gamma, DTASpeed, DTAP, TDAvgSpeedDiff);

                ///////////////////////////////////////  noon
                peak_no = 1;

                assign_period_start_time_in_hour = 10;
                assign_period_end_time_in_hour = 14;
                starting_time_in_hour = 12;
                ending_time_in_hour = 13;
                t2 = 12.5;  //12 noon

                obs_P_in_hour = g_TMC_vector[tmc_index].scan_congestion_duration(peak_no, starting_time_in_hour,
                    ending_time_in_hour, FD_vc, obs_t0_in_hour, obs_t3_in_hour, pLink, assignment_D, assignment_VMT, assignment_VHT, assignment_VDT, assignment_VCDT, congestion_D,
                    congestion_dc_ratio, mean_congestion_mu, mean_congestion_speed, highest_speed, mean_speed, t2_speed, t2_queue, gamma);

                V = assignment_D;
                vcd = pLink->vc;
                vct = pLink->vc;

                g_TMC_vector[tmc_index].perform_estimation(false, pLink, peak_no, assign_period_start_time_in_hour, assign_period_end_time_in_hour, t2, V, laneCapacity, vf, vcd, vct, DTASpeed, DTAP, TDAvgSpeedDiff);
                STA_MD_MAE = fabs(g_link_vector[i].VDF_STA_speed[1] - mean_speed);
                STA_MD_APE = STA_MD_MAE / max(1, mean_speed);

                QHF = 3;
                if (congestion_D > 1 && obs_P_in_hour >= 1)
                    QHF = V / congestion_D;

                fprintf(p_file_tmc_link, "%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,",
                    FD_vc,
                    QHF, g_TMC_vector[tmc_index].est_QDF_n[peak_no], g_TMC_vector[tmc_index].est_QDF_s[peak_no],
                    obs_t0_in_hour, obs_t3_in_hour, obs_P_in_hour, assignment_D, assignment_VMT, assignment_VHT, assignment_VDT, assignment_VCDT,
                    congestion_D, congestion_dc_ratio, mean_congestion_mu, mean_congestion_speed,
                    highest_speed, mean_speed,
                    t2_speed, t2_queue, gamma, DTASpeed, DTAP, TDAvgSpeedDiff);

                ///////////////////////////////////////  PM
                peak_no = 2;

                assign_period_start_time_in_hour = 14;
                assign_period_end_time_in_hour = 20;
                starting_time_in_hour = 15;
                ending_time_in_hour = 19;
                t2 = 17;  //5pm

                obs_P_in_hour = g_TMC_vector[tmc_index].scan_congestion_duration(peak_no, starting_time_in_hour,
                    ending_time_in_hour, FD_vc, obs_t0_in_hour, obs_t3_in_hour, pLink, assignment_D, assignment_VMT, assignment_VHT, assignment_VDT, assignment_VCDT, congestion_D,
                    congestion_dc_ratio, mean_congestion_mu, mean_congestion_speed, highest_speed, mean_speed, t2_speed, t2_queue, gamma);

                V = assignment_D;
                vcd = pLink->vc;
                vct = pLink->vc;

                g_TMC_vector[tmc_index].perform_estimation(false, pLink, peak_no, assign_period_start_time_in_hour, assign_period_end_time_in_hour, t2, V, laneCapacity, vf, vcd, vct, DTASpeed, DTAP, TDAvgSpeedDiff);
                STA_PM_MAE = fabs(g_link_vector[i].VDF_STA_speed[3] - mean_speed);
                STA_PM_APE = STA_PM_MAE / max(1, mean_speed);

                QHF = 3;
                if (congestion_D > 1 && obs_P_in_hour >= 1)
                    QHF = max(1, V / congestion_D);

                fprintf(p_file_tmc_link, "%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,",
                    FD_vc,
                    QHF, g_TMC_vector[tmc_index].est_QDF_n[peak_no], g_TMC_vector[tmc_index].est_QDF_s[peak_no],
                    obs_t0_in_hour, obs_t3_in_hour, obs_P_in_hour, assignment_D, assignment_VMT, assignment_VHT, assignment_VDT, assignment_VCDT,
                    congestion_D, congestion_dc_ratio, mean_congestion_mu, mean_congestion_speed,
                    highest_speed, mean_speed,
                    t2_speed, t2_queue, gamma, DTASpeed, DTAP, TDAvgSpeedDiff);

            }
            else
            {
                // no data AM

                fprintf(p_file_tmc_link, "%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,",
                    FD_vc, -1, -1, -1, obs_t0_in_hour, obs_t3_in_hour, obs_P_in_hour, assignment_D, assignment_VMT, assignment_VHT, assignment_VDT, assignment_VCDT,
                    congestion_D, congestion_dc_ratio, mean_congestion_mu, mean_congestion_speed,
                    highest_speed, mean_speed,
                    t2_speed, t2_queue, gamma, DTASpeed, DTAP, TDAvgSpeedDiff);
                // no data MD
                fprintf(p_file_tmc_link, "%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,",
                    FD_vc, -1, -1, -1, obs_t0_in_hour, obs_t3_in_hour, obs_P_in_hour, assignment_D, assignment_VMT, assignment_VHT, assignment_VDT, assignment_VCDT,
                    congestion_D, congestion_dc_ratio, mean_congestion_mu, mean_congestion_speed,
                    highest_speed, mean_speed,
                    t2_speed, t2_queue, gamma, DTASpeed, DTAP, TDAvgSpeedDiff);

                // no data PM
                fprintf(p_file_tmc_link, "%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,",
                    FD_vc, -1, -1, -1, obs_t0_in_hour, obs_t3_in_hour, obs_P_in_hour, assignment_D, assignment_VMT, assignment_VHT, assignment_VDT, assignment_VCDT,
                    congestion_D, congestion_dc_ratio, mean_congestion_mu, mean_congestion_speed,
                    highest_speed, mean_speed,
                    t2_speed, t2_queue, gamma, DTASpeed, DTAP, TDAvgSpeedDiff);
            }
            // vc, t0, t3, P, D, mu,

            fprintf(p_file_tmc_link, "\"%s\",", g_link_vector[i].geometry.c_str());

            fprintf(p_file_tmc_link, "\"LINESTRING (");

            fprintf(p_file_tmc_link, "%f %f,", g_link_vector[i].TMC_from.x, g_link_vector[i].TMC_from.y);
            fprintf(p_file_tmc_link, "%f %f,", g_link_vector[i].TMC_to.x, g_link_vector[i].TMC_to.y);
            fprintf(p_file_tmc_link, ")\"");
            fprintf(p_file_tmc_link, ",");

            /* if (bReadingDataReady)*/
            {
                if (assignment.m_TMClink_map.find(g_link_vector[i].tmc_code) != assignment.m_TMClink_map.end())
                {
                    int tmc_index = assignment.m_TMClink_map[g_link_vector[i].tmc_code];

                    for (int tau = 1; tau <= 4; tau++)
                    {
                        fprintf(p_file_tmc_link, "%f,%f,%f,",
                            g_link_vector[i].VDF_STA_speed[tau], g_link_vector[i].VDF_STA_VOC[tau], g_link_vector[i].VDF_STA_volume[tau]);
                    }

                    double ObsSpeed[25], EstSpeed[25], EstSpeedDiff[25];
                    double MAE_total = 0;
                    double MAPE_total = 0;
                    double RMSE_total = 0;
                    int count_total = 0;

                    for (int t = 6 * 60; t < 20 * 60; t += 60)
                    {
                        int hour = t / 60;
                        ObsSpeed[hour] = g_TMC_vector[tmc_index].get_avg_hourly_speed(t);
                        fprintf(p_file_tmc_link, "%.1f,", ObsSpeed[hour]);
                    }
                    for (int t = 6 * 60; t < 20 * 60; t += 60)
                    {
                        int hour = t / 60;
                        EstSpeed[hour] = g_TMC_vector[tmc_index].get_est_hourly_speed(t);
                        if (EstSpeed[hour] > 1)  // valid
                        {
                            EstSpeedDiff[hour] = fabs(EstSpeed[hour] - ObsSpeed[hour]);
                            MAE_total += fabs(EstSpeedDiff[hour]);
                            MAPE_total += fabs(EstSpeedDiff[hour]) / max(1, ObsSpeed[hour]);
                            RMSE_total += EstSpeedDiff[hour] * EstSpeedDiff[hour];
                            count_total += 1;
                        }
                        else
                            EstSpeedDiff[hour] = -1;
                    }
                    for (int t = 6 * 60; t < 20 * 60; t += 60)
                    {
                        int hour = t / 60;
                        EstSpeed[hour] = g_TMC_vector[tmc_index].get_est_hourly_speed(t);

                        fprintf(p_file_tmc_link, "%.1f,", EstSpeed[hour]);
                    }
                    for (int t = 6 * 60; t < 20 * 60; t += 60)
                    {
                        int hour = t / 60;
                        fprintf(p_file_tmc_link, "%.1f,", EstSpeedDiff[hour]);
                    }

                    double MSE_total = RMSE_total / max(1, count_total);

                    fprintf(p_file_tmc_link, "%.2f,%.2f,%.2f,", MAE_total / max(1, count_total), MAPE_total / max(1, count_total) * 100, pow(MSE_total, 0.5));

                    /// <summary>
                    ///
                    /// </summary>
                    /// <param name="bReadingDataReady"></param>
                    for (int t = 6 * 60; t < 20 * 60; t += 60)
                    {
                        int hour = t / 60;
                        ObsSpeed[hour] = g_TMC_vector[tmc_index].get_avg_hourly_speed(t);

                        float speed_ratio = ObsSpeed[hour] / max(1, g_link_vector[i].TMC_highest_speed);

                        if (speed_ratio >= 1)
                            speed_ratio = 1;
                        fprintf(p_file_tmc_link, "%.3f,", speed_ratio);
                    }
                    for (int t = 6 * 60; t < 20 * 60; t += 60)
                    {
                        int hour = t / 60;
                        EstSpeed[hour] = g_TMC_vector[tmc_index].get_est_hourly_speed(t);
                        float speed_ratio = EstSpeed[hour] / max(1, g_link_vector[i].TMC_highest_speed);

                        if (speed_ratio >= 1)
                            speed_ratio = 1;
                        fprintf(p_file_tmc_link, "%.3f,", speed_ratio);
                    }

                    double STA_RMSE = pow((STA_AM_MAE * STA_AM_MAE + STA_MD_MAE * STA_MD_MAE + STA_PM_MAE * STA_PM_MAE) / 3, 0.5);
                    fprintf(p_file_tmc_link, "%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,", STA_AM_MAE, STA_MD_MAE, STA_PM_MAE, (STA_AM_MAE + STA_MD_MAE + STA_PM_MAE) / 3, (STA_AM_APE + STA_MD_APE + STA_PM_APE) / 3 * 100, STA_RMSE);

                    for (int t = 6 * 60; t < 20 * 60; t += 5)
                    {

                        fprintf(p_file_tmc_link, "%.1f,", g_TMC_vector[tmc_index].get_avg_speed(t));
                    }

                    for (int t = 0 * 60; t < 24 * 60; t += 15)
                    {

                        fprintf(p_file_tmc_link, "%.1f,", g_TMC_vector[tmc_index].get_avg_speed_15min(t));
                    }

                    for (int t = 6 * 60; t < 20 * 60; t += 15)
                    {
                        float speed_ratio = g_TMC_vector[tmc_index].get_avg_speed_15min(t) / max(1, g_link_vector[i].TMC_highest_speed);

                        if (speed_ratio >= 1)
                            speed_ratio = 1;

                        fprintf(p_file_tmc_link, "%.2f,", speed_ratio);
                    }

                    for (int t = 6 * 60; t < 20 * 60; t += 15)
                    {
                        float speed = g_TMC_vector[tmc_index].get_avg_speed(t);
                        double volume = g_link_vector[i].get_volume_from_speed(speed, g_link_vector[i].TMC_highest_speed);

                        fprintf(p_file_tmc_link, "%.2f,", volume);
                    }

                    for (int t = 6 * 60; t < 20 * 60; t += 15)
                    {
                        float speed = g_TMC_vector[tmc_index].get_avg_speed_15min(t);
                        double volume = g_link_vector[i].get_volume_from_speed(speed, g_link_vector[i].TMC_highest_speed);

                        fprintf(p_file_tmc_link, "%.2f,", volume);
                    }

                    // estimation 
                    for (int t = 6 * 60; t < 20 * 60; t += 15)
                    {

                        fprintf(p_file_tmc_link, "%.1f,", g_TMC_vector[tmc_index].get_est_speed(t));
                    }

                    //    for (int t = 6 * 60; t < 20 * 60; t += 5)
                    //    {
                    //        float est_speed = g_TMC_vector[tmc_index].get_est_speed(t);
                    //        float speed_ratio = 1;

                    //        if(est_speed >1)  // valid result
                    //        speed_ratio = g_TMC_vector[tmc_index].get_est_speed(t) / max(1, g_link_vector[i].TMC_highest_speed);

                    //        if (speed_ratio >= 1)
                    //            speed_ratio = 1;


                    //        fprintf(p_file_tmc_link, "%.2f,", speed_ratio);
                    //    }


                }
            }
            fprintf(p_file_tmc_link, "\n");

            if (assignment.m_TMClink_map.find(g_link_vector[i].tmc_code) != assignment.m_TMClink_map.end())
            {
                int tmc_index = assignment.m_TMClink_map[g_link_vector[i].tmc_code];
                // output link QVDF parameters 
                fprintf(p_file_tmc_vdf, "%s,%d,%d,%s,%.2f,%.2f,%.2f,",
                    g_link_vector[i].link_id.c_str(),
                    g_node_vector[g_link_vector[i].from_node_seq_no].node_id,
                    g_node_vector[g_link_vector[i].to_node_seq_no].node_id,
                    g_link_vector[i].link_type_code.c_str(),
                    g_link_vector[i].free_speed,
                    g_link_vector[i].lane_capacity,
                    g_link_vector[i].kc
                );

                //            fprintf(p_file_tmc_vdf, "QVDF_qhf%d,QVDF_n%d,QVDF_s%d,QVDF_v%d,QVDF_d%d,QVDF_dc%d,QVDF_p%d,QVDF_mu%d,QVDF_vcd%d,",

                for (int k = 0; k < 3; k++)
                {

                    fprintf(p_file_tmc_vdf, "%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,",
                        g_TMC_vector[tmc_index].est_PLF[k],
                        g_link_vector[i].s3_m,
                        g_TMC_vector[tmc_index].est_QDF_n[k],
                        g_TMC_vector[tmc_index].est_QDF_s[k],
                        g_TMC_vector[tmc_index].est_QDF_beta[k],
                        g_TMC_vector[tmc_index].est_gamma[k],
                        g_TMC_vector[tmc_index].est_V[k],
                        g_TMC_vector[tmc_index].est_D[k],
                        g_TMC_vector[tmc_index].est_DCRatio[k],
                        g_TMC_vector[tmc_index].est_P[k],
                        g_TMC_vector[tmc_index].est_mu[k],
                        g_TMC_vector[tmc_index].est_vcd[k]);
               }

            }

            fprintf(p_file_tmc_vdf, "\n");
        }



        

        fclose(p_file_tmc_link);
        fclose(p_file_tmc_vdf);
    }
    else
    {
        dtalog.output() << "Error: File link_cbi_summary.csv cannot be opened.\n It might be currently used and locked by EXCEL." << endl;
        g_program_stop();

    }

}
void g_output_tmc_scenario_files()
{

    CCSVParser parser_scenario;

    if (parser_scenario.OpenCSVFile("link_scenario.csv", true))
    {
        while (parser_scenario.ReadRecord())  // if this line contains [] mark, then we will also read field headers.
        {
            //string tmc;

            //for (int i = 0; i < g_link_vector.size(); i++)
            //{
            //    if (g_link_vector[i].link_type >= 0 && g_link_vector[i].tmc_code.size() > 0)
            //    {
            //        if (g_link_vector[i].tmc_code.compare(tmc)!=0)
            //        {
            //            g_link_vector[i].Scenario_evaluation_flag = true;

            //            parser_scenario.GetValueByFieldName("QDF_voc_ratio1", g_link_vector[i].Scenario_STA_VOC_Ratio[1], false);
            //            parser_scenario.GetValueByFieldName("QDF_voc_ratio2", g_link_vector[i].Scenario_STA_VOC_Ratio[2], false);
            //            parser_scenario.GetValueByFieldName("QDF_voc_ratio3", g_link_vector[i].Scenario_STA_VOC_Ratio[3], false);
            //            parser_scenario.GetValueByFieldName("QDF_voc_ratio4", g_link_vector[i].Scenario_STA_VOC_Ratio[4], false);



            //        }
            //    }
            //}

            int from_node_id;
            if (!parser_scenario.GetValueByFieldName("from_node_id", from_node_id))
                continue;

            int to_node_id;
            if (!parser_scenario.GetValueByFieldName("to_node_id", to_node_id))
                continue;


            if (assignment.g_node_id_to_seq_no_map.find(from_node_id) == assignment.g_node_id_to_seq_no_map.end())
            {
                dtalog.output() << "Error: from_node_id " << from_node_id << " in file scenario.csv is not defined in node.csv." << endl;
                continue; //has not been defined
            }

            if (assignment.g_node_id_to_seq_no_map.find(to_node_id) == assignment.g_node_id_to_seq_no_map.end())
            {
                dtalog.output() << "Error: to_node_id " << to_node_id << " in file scenario.csv is not defined in node.csv." << endl;
                continue; //has not been defined
            }


            int internal_from_node_seq_no = assignment.g_node_id_to_seq_no_map[from_node_id];  // map external node number to internal node seq no.
            int internal_to_node_seq_no = assignment.g_node_id_to_seq_no_map[to_node_id];

            if (g_node_vector[internal_from_node_seq_no].m_to_node_2_link_seq_no_map.find(internal_to_node_seq_no) != g_node_vector[internal_from_node_seq_no].m_to_node_2_link_seq_no_map.end())
            {
                int link_seq_no = g_node_vector[internal_from_node_seq_no].m_to_node_2_link_seq_no_map[internal_to_node_seq_no];
                int i = link_seq_no;
                g_link_vector[i].Scenario_evaluation_flag = true;

                parser_scenario.GetValueByFieldName("QDF_voc_ratio1", g_link_vector[i].Scenario_STA_VOC_Ratio[1], false);
                parser_scenario.GetValueByFieldName("QDF_voc_ratio2", g_link_vector[i].Scenario_STA_VOC_Ratio[2], false);
                parser_scenario.GetValueByFieldName("QDF_voc_ratio3", g_link_vector[i].Scenario_STA_VOC_Ratio[3], false);
                parser_scenario.GetValueByFieldName("QDF_voc_ratio4", g_link_vector[i].Scenario_STA_VOC_Ratio[4], false);


            }
        }

    }

    FILE* p_file_tmc_link = fopen("link_cbsc.csv", "w");

    if (p_file_tmc_link != NULL)
    {
        fprintf(p_file_tmc_link, "link_id,tmc,tmc_corridor_name,tmc_corridor_id,tmc_road_order,tmc_road_sequence,tmc_road,tmc_direction,tmc_intersection,link_no,from_node_id,to_node_id,link_type,");
        fprintf(p_file_tmc_link, "link_type_code,FT,AT,lanes,free_speed,capacity, kc, ");
        fprintf(p_file_tmc_link, "AM_vc, AM_t0, AM_t3, AM_P, AM_Assign_V, AM_D,AM_DC_ratio_base,AM_DC_ratio_future, AM_mu, AM_vu, AM_vf_reference, AM_v_mean, AM_t2_v, AM_t2_queue, AM_gamma, AM_DTASpeed1, AM_DTAP1, AM_DTATDSpdDiff,");
        fprintf(p_file_tmc_link, "PM_vc, PM_t0,PM_t3,PM_P,PM_Assign_V,PM_D,PM_DC_ratio_base,PM_DC_ratio_future,PM_mu,PM_vu,PM_vf_reference,PM_v_mean, PM_t2_v, PM_t2_queue, PM_gamma, PM_DTASpeed1, PM_DTAP1, PM_DTATDSpdDiff,");
        fprintf(p_file_tmc_link, ",geometry,tmc_geometry,reading_road_name,QDM_voc_ratio1,QDM_voc_ratio2,QDM_voc_ratio3,");

        fprintf(p_file_tmc_link, "|,");

        for (int t = 6 * 60; t < 20 * 60; t += 60)
        {
            int hour = t / 60;
            int minute = t - hour * 60;

            fprintf(p_file_tmc_link, "vh%02d,", hour, minute);
        }
        fprintf(p_file_tmc_link, "|,");
        for (int t = 6 * 60; t < 20 * 60; t += 60)
        {
            int hour = t / 60;

            fprintf(p_file_tmc_link, "pvh%02d,", hour);
        }
        fprintf(p_file_tmc_link, "|,");

        for (int t = 6 * 60; t < 20 * 60; t += 60)
        {
            int hour = t / 60;

            fprintf(p_file_tmc_link, "pvh%02ddiff,", hour);
        }
        fprintf(p_file_tmc_link, "|,");

        fprintf(p_file_tmc_link, "pvhAE,pvhAPE,");

        fprintf(p_file_tmc_link, "|,");
        // travel time 
        for (int t = 6 * 60; t < 20 * 60; t += 60)
        {
            int hour = t / 60;

            fprintf(p_file_tmc_link, "tth%02d,", hour);
        }
        fprintf(p_file_tmc_link, "|,");
        for (int t = 6 * 60; t < 20 * 60; t += 60)
        {
            int hour = t / 60;

            fprintf(p_file_tmc_link, "ptth%02d,", hour);
        }
        fprintf(p_file_tmc_link, "|,");

        // travel time 
        for (int t = 6 * 60; t < 20 * 60; t += 60)
        {
            int hour = t / 60;

            fprintf(p_file_tmc_link, "vtth%02d,", hour);
        }
        fprintf(p_file_tmc_link, "|,");
        for (int t = 6 * 60; t < 20 * 60; t += 60)
        {
            int hour = t / 60;

            fprintf(p_file_tmc_link, "ptth%02d,", hour);
        }
        fprintf(p_file_tmc_link, "|,");

        // 5 min prediction 

        for (int t = 6 * 60; t < 20 * 60; t += 5)
        {
            int hour = t / 60;
            int minute = t - hour * 60;

            fprintf(p_file_tmc_link, "tt%02d:%02d,", hour, minute);
        }

        fprintf(p_file_tmc_link, "\n");

        std::map<_int64, int> TMC_long_id_mapping;  // this is used to mark if this cell_id has been identified or not

        //// sort data records
        for (int i = 0; i < g_link_vector.size(); i++)
        {
            if (g_link_vector[i].link_type >= 0 && g_link_vector[i].tmc_code.size() > 0)
            {
                _int64 TMC_long_key = (g_link_vector[i].tmc_corridor_id * 10000 + g_link_vector[i].tmc_road_sequence) * 1000000 + g_link_vector[i].link_seq_no;
                TMC_long_id_mapping[TMC_long_key] = g_link_vector[i].link_seq_no;
            }
        }


        std::map<_int64, int>::iterator it;

        for (it = TMC_long_id_mapping.begin(); it != TMC_long_id_mapping.end(); ++it)
        {
            int i = it->second;

            if (g_link_vector[i].link_type >= 0 && g_link_vector[i].tmc_code.size() > 0)
            {

                if (assignment.m_TMClink_map.find(g_link_vector[i].tmc_code) != assignment.m_TMClink_map.end())
                {
                    int tmc_index = assignment.m_TMClink_map[g_link_vector[i].tmc_code];
                    if (g_TMC_vector[tmc_index].b_with_sensor_speed_data == false)  // with reading speed data only
                        continue;
                }



                float free_speed = g_link_vector[i].free_speed;
                g_link_vector[i].update_kc(g_link_vector[i].free_speed);

                fprintf(p_file_tmc_link, "%s,%s,%s,%d,%d,%d,%s,%s,%s,%d,%d,%d,%d,%s,%d,%d,%d,%f,%f,%f,",
                    g_link_vector[i].link_id.c_str(),
                    g_link_vector[i].tmc_code.c_str(),
                    g_link_vector[i].tmc_corridor_name.c_str(),
                    g_link_vector[i].tmc_corridor_id,
                    g_link_vector[i].tmc_road_order,
                    g_link_vector[i].tmc_road_sequence,
                    g_link_vector[i].tmc_road.c_str(),
                    g_link_vector[i].tmc_direction.c_str(),
                    g_link_vector[i].tmc_intersection.c_str(),
                    g_link_vector[i].link_seq_no,
                    g_node_vector[g_link_vector[i].from_node_seq_no].node_id,
                    g_node_vector[g_link_vector[i].to_node_seq_no].node_id,
                    g_link_vector[i].link_type,
                    g_link_vector[i].link_type_code.c_str(),
                    g_link_vector[i].FT,
                    g_link_vector[i].AT,
                    g_link_vector[i].number_of_lanes,
                    g_link_vector[i].free_speed,
                    g_link_vector[i].lane_capacity,
                    g_link_vector[i].kc

                );

                float FD_vc = 0;
                float obs_t0_in_hour = -1;
                float obs_t3_in_hour = -1;
                float mean_congestion_speed = 0;
                float congestion_dc_ratio = 0;
                float obs_P_in_hour = 0;
                float mean_congestion_mu = 0;
                float highest_speed = 0;
                float mean_speed = 0;
                float congestion_D = 0;
                float assignment_D = 0;
                float assignment_VMT = 0;
                float assignment_VHT = 0;
                float assignment_VDT = 0;
                float assignment_VCDT = 0;
                float t2_speed = 0;
                float t2_queue = 0;
                float gamma = 0;
                float assign_period_start_time_in_hour = 6;
                float assign_period_end_time_in_hour = 11;
                int starting_time_in_hour = 6;
                int ending_time_in_hour = 9;
                float t2 = 7.5;
                float DTASpeed = -1;
                float DTAP = -1;
                float TDAvgSpeedDiff = -1;

                // P D analysis
                if (assignment.m_TMClink_map.find(g_link_vector[i].tmc_code) != assignment.m_TMClink_map.end())
                {
                    int tmc_index = assignment.m_TMClink_map[g_link_vector[i].tmc_code];

                    int peak_no = 0;  //AM

                    assign_period_start_time_in_hour = 6;
                    assign_period_end_time_in_hour = 10;
                    starting_time_in_hour = 6;
                    ending_time_in_hour = 9;

                    congestion_dc_ratio = 0;

                    CLink* pLink = &(g_link_vector[i]);

                    int from_node_id = g_node_vector[g_link_vector[i].from_node_seq_no].node_id;
                    int to_node_id = g_node_vector[g_link_vector[i].to_node_seq_no].node_id;

                    if (from_node_id == 19229 && to_node_id == 19233)
                    {
                        int idebug = 1;
                    }

                    obs_P_in_hour = g_TMC_vector[tmc_index].scan_congestion_duration(peak_no, starting_time_in_hour,
                        ending_time_in_hour, FD_vc, obs_t0_in_hour, obs_t3_in_hour, pLink, assignment_D, assignment_VMT, assignment_VHT, assignment_VDT, assignment_VCDT, congestion_D,
                        congestion_dc_ratio, mean_congestion_mu, mean_congestion_speed, highest_speed, mean_speed, t2_speed, t2_queue, gamma);


                    float V = assignment_D;  // apply demand flow change: AM
                    float laneCapacity = pLink->lane_capacity;
                    float vf = highest_speed;
                    float vcd = pLink->vc;
                    float vct = pLink->vc;
                    t2 = 7.5;

                    if (pLink->Scenario_evaluation_flag && g_link_vector[i].tmc_code.compare("115 + 04402") != 0)
                    {
                        int i_debug = 1;
                    }


                    if (pLink->Scenario_STA_VOC_Ratio[1] > 1.0001)
                    {
                        int i_debug = 1;
                    }

                    float AM_DC_ratio_base, AM_DC_ratio_future;
                    g_TMC_vector[tmc_index].perform_estimation(false, pLink, peak_no, assign_period_start_time_in_hour, assign_period_end_time_in_hour, t2, V,
                        laneCapacity, vf, vcd, vct, DTASpeed, DTAP, TDAvgSpeedDiff);

                    AM_DC_ratio_base = g_TMC_vector[tmc_index].est_DCRatio[peak_no];

                    V = assignment_D * pLink->Scenario_STA_VOC_Ratio[1];  // apply demand flow change: AM

                    g_TMC_vector[tmc_index].perform_estimation(true, pLink, peak_no, assign_period_start_time_in_hour, assign_period_end_time_in_hour, t2, V,
                        laneCapacity, vf, vcd, vct, DTASpeed, DTAP, TDAvgSpeedDiff);

                    AM_DC_ratio_future = g_TMC_vector[tmc_index].est_DCRatio[peak_no];

                    fprintf(p_file_tmc_link, "%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,",
                        FD_vc, obs_t0_in_hour, obs_t3_in_hour, obs_P_in_hour, assignment_D,
                        congestion_D, AM_DC_ratio_base, AM_DC_ratio_future, g_TMC_vector[tmc_index].est_mu[peak_no], mean_congestion_speed,
                        highest_speed, mean_speed,
                        g_TMC_vector[tmc_index].est_vt2[peak_no], t2_queue, gamma, DTASpeed, DTAP, TDAvgSpeedDiff);

                    ///////////////////////////////////////  noon
                    peak_no = 1;

                    assign_period_start_time_in_hour = 10;
                    assign_period_end_time_in_hour = 14;
                    starting_time_in_hour = 12;
                    ending_time_in_hour = 13;
                    t2 = 12;  //12 noon


                    obs_P_in_hour = g_TMC_vector[tmc_index].scan_congestion_duration(peak_no, starting_time_in_hour,
                        ending_time_in_hour, FD_vc, obs_t0_in_hour, obs_t3_in_hour, pLink, assignment_D, assignment_VMT, assignment_VHT, assignment_VDT, assignment_VCDT, congestion_D,
                        congestion_dc_ratio, mean_congestion_mu, mean_congestion_speed, highest_speed, mean_speed, t2_speed, t2_queue, gamma);

                    V = assignment_D;
                    vcd = pLink->vc;
                    vct = pLink->vc;

                    g_TMC_vector[tmc_index].perform_estimation(false, pLink, peak_no, assign_period_start_time_in_hour, assign_period_end_time_in_hour, t2, V, laneCapacity, vf, vcd, vct, DTASpeed, DTAP, TDAvgSpeedDiff);
                    V = assignment_D * pLink->Scenario_STA_VOC_Ratio[2];  // apply demand flow change noon
                    g_TMC_vector[tmc_index].perform_estimation(true, pLink, peak_no, assign_period_start_time_in_hour, assign_period_end_time_in_hour, t2, V, laneCapacity, vf, vcd, vct, DTASpeed, DTAP, TDAvgSpeedDiff);

                    ///////////////////////////////////////  PM
                    peak_no = 2;

                    assign_period_start_time_in_hour = 14;
                    assign_period_end_time_in_hour = 20;
                    starting_time_in_hour = 15;
                    ending_time_in_hour = 19;
                    t2 = 17;  //5pm

                    obs_P_in_hour = g_TMC_vector[tmc_index].scan_congestion_duration(peak_no, starting_time_in_hour,
                        ending_time_in_hour, FD_vc, obs_t0_in_hour, obs_t3_in_hour, pLink, assignment_D, assignment_VMT, assignment_VHT, assignment_VDT, assignment_VCDT, congestion_D,
                        congestion_dc_ratio, mean_congestion_mu, mean_congestion_speed, highest_speed, mean_speed, t2_speed, t2_queue, gamma);

                    V = assignment_D;
                    vcd = pLink->vc;
                    vct = pLink->vc;

                    float PM_DC_ratio_base, PM_DC_ratio_future;

                    g_TMC_vector[tmc_index].perform_estimation(false, pLink, peak_no, assign_period_start_time_in_hour, assign_period_end_time_in_hour, t2, V, laneCapacity, vf, vcd, vct, DTASpeed, DTAP, TDAvgSpeedDiff);

                    PM_DC_ratio_base = g_TMC_vector[tmc_index].est_DCRatio[peak_no];
                    V = assignment_D * pLink->Scenario_STA_VOC_Ratio[3];  // apply demand flow change noon
                    g_TMC_vector[tmc_index].perform_estimation(true, pLink, peak_no, assign_period_start_time_in_hour, assign_period_end_time_in_hour, t2, V, laneCapacity, vf, vcd, vct, DTASpeed, DTAP, TDAvgSpeedDiff);
                    PM_DC_ratio_future = g_TMC_vector[tmc_index].est_DCRatio[peak_no];
                    fprintf(p_file_tmc_link, "%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,",
                        FD_vc, obs_t0_in_hour, obs_t3_in_hour, obs_P_in_hour, assignment_D,
                        congestion_D, PM_DC_ratio_base, PM_DC_ratio_future, g_TMC_vector[tmc_index].est_mu[peak_no], mean_congestion_speed,
                        highest_speed, mean_speed,
                        g_TMC_vector[tmc_index].est_vt2[peak_no], t2_queue, gamma, DTASpeed, DTAP, TDAvgSpeedDiff);
                }
                else
                {
                    // no data

                    fprintf(p_file_tmc_link, "%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,",
                        FD_vc, obs_t0_in_hour, obs_t3_in_hour, obs_P_in_hour, assignment_D,
                        congestion_D, congestion_dc_ratio, mean_congestion_mu, mean_congestion_speed,
                        highest_speed, mean_speed,
                        t2_speed, t2_queue, gamma, DTASpeed, DTAP, TDAvgSpeedDiff);
                }
                // vc, t0, t3, P, D, mu,

                fprintf(p_file_tmc_link, "\"%s\",", g_link_vector[i].geometry.c_str());

                fprintf(p_file_tmc_link, "\"LINESTRING (");

                fprintf(p_file_tmc_link, "%f %f,", g_link_vector[i].TMC_from.x, g_link_vector[i].TMC_from.y);
                fprintf(p_file_tmc_link, "%f %f,", g_link_vector[i].TMC_to.x, g_link_vector[i].TMC_to.y);
                fprintf(p_file_tmc_link, ")\"");
                fprintf(p_file_tmc_link, ",,");
                fprintf(p_file_tmc_link, "%f,%f,%f,", g_link_vector[i].Scenario_STA_VOC_Ratio[1], g_link_vector[i].Scenario_STA_VOC_Ratio[2], g_link_vector[i].Scenario_STA_VOC_Ratio[3]);


                if (assignment.m_TMClink_map.find(g_link_vector[i].tmc_code) != assignment.m_TMClink_map.end())
                {
                    int tmc_index = assignment.m_TMClink_map[g_link_vector[i].tmc_code];

                    double PredSpeed[25], EstSpeed[25], PredSpeedDiff[25];
                    double PredTravelTime[25], EstTravelTime[25], PredTravelTimeDiff[25];
                    double MAE_total = 0;
                    double TravelTime_total = 0;
                    double MAPE_total = 0;
                    int count_total = 0;

                    for (int t = 6 * 60; t < 20 * 60; t += 60)
                    {
                        int hour = t / 60;
                        EstSpeed[hour] = g_TMC_vector[tmc_index].get_est_hourly_speed(t);
                        fprintf(p_file_tmc_link, "%.1f,", EstSpeed[hour]);
                    }
                    fprintf(p_file_tmc_link, ",");
                    for (int t = 6 * 60; t < 20 * 60; t += 60)
                    {
                        int hour = t / 60;
                        PredSpeed[hour] = g_TMC_vector[tmc_index].get_pred_hourly_speed(t);
                        PredSpeedDiff[hour] = PredSpeed[hour] - EstSpeed[hour];
                        MAE_total += PredSpeedDiff[hour];
                        MAPE_total += fabs(PredSpeedDiff[hour]) / max(1, PredSpeed[hour]);

                        count_total += 1;

                        fprintf(p_file_tmc_link, "%.1f,", PredSpeed[hour]);
                    }
                    fprintf(p_file_tmc_link, ",");
                    for (int t = 6 * 60; t < 20 * 60; t += 60)
                    {
                        int hour = t / 60;
                        fprintf(p_file_tmc_link, "%.1f,", PredSpeedDiff[hour]);
                    }
                    fprintf(p_file_tmc_link, ",");

                    fprintf(p_file_tmc_link, "%.2f,%.2f,", MAE_total / max(1, count_total), MAPE_total / max(1, count_total) * 100);
                    fprintf(p_file_tmc_link, ",");

                    // travel time, to be computed for corridor level travel time
                    for (int t = 6 * 60; t < 20 * 60; t += 60)
                    {
                        int hour = t / 60;
                        EstSpeed[hour] = g_TMC_vector[tmc_index].get_est_hourly_speed(t);
                        EstTravelTime[hour] = g_link_vector[i].link_distance_in_km / max(1, EstSpeed[hour]) * 60;

                        fprintf(p_file_tmc_link, "%.2f,", EstTravelTime[hour]);
                    }
                    fprintf(p_file_tmc_link, ",");
                    for (int t = 6 * 60; t < 20 * 60; t += 60)
                    {
                        int hour = t / 60;
                        PredSpeed[hour] = g_TMC_vector[tmc_index].get_pred_hourly_speed(t);
                        PredTravelTime[hour] = g_link_vector[i].link_distance_in_km / max(1, PredSpeed[hour]) * 60;
                        count_total += 1;

                        fprintf(p_file_tmc_link, "%.2f,", PredTravelTime[hour]);
                    }

                    fprintf(p_file_tmc_link, ",");

                    // 
                    for (int t = 6 * 60; t < 20 * 60; t += 5)
                    {
                        fprintf(p_file_tmc_link, "%.1f,", g_link_vector[i].link_distance_in_km / max(1, g_TMC_vector[tmc_index].pred_speed[t]) * 60);
                    }
                }
                fprintf(p_file_tmc_link, "\n");

            }


        }

        fclose(p_file_tmc_link);
    }
    else
    {
        dtalog.output() << "Error: File TMC_link_scenario.csv cannot be opened.\n It might be currently used and locked by EXCEL." << endl;
        g_program_stop();

    }

}
