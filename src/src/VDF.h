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
//using std::istringstream;

class CPeriod_VDF
{
public:
    CPeriod_VDF() : n{ 1.28 }, s{ 0.859 }, vc{ 45 }, m{ 0.5 }, peak_load_factor { 1 }, VOC{ 0 }, gamma{ 3.47f }, mu{ 1000 }, PHF{ 3 },
        alpha{ 0.15f }, beta{ 4 }, rho{ 1 }, preload{ 0 }, penalty{ 0 }, LR_price{ 0 }, LR_RT_price{ 0 }, marginal_base{ 1 },
        starting_time_slot_no{ 0 }, ending_time_slot_no{ 0 },
        cycle_link_distance_in_km{ -1 }, red_time{ 0 }, effective_green_time{ 0 }, t0{ 0 }, t3{ 0 }, start_green_time{ -1 }, end_green_time{ -1 }, sav{ 0 }, time_period_in_hour{ 1 }
    {
        for (int at = 0; at < MAX_AGNETTYPES; at++)
        {
            toll[at] = 0;
            pce[at] = 1;
        }

    }

    ~CPeriod_VDF()
    {
    }

    //float PerformSignalVDF(float hourly_per_lane_volume, float red, float cycle_link_distance_in_km)
    //{
    //    float lambda = hourly_per_lane_volume;
    //    float mu = _default_saturation_flow_rate; //default saturation flow ratesa
    //    float s_bar = 1.0 / 60.0 * red * red / (2 * cycle_link_distance_in_km); // 60.0 is used to convert sec to min
    //    float uniform_delay = s_bar / max(1 - lambda / mu, 0.1f);

    //    return uniform_delay;
    //}

    double calculate_travel_time_based_on_VDF(double volume, int VDF_type_no, int lanes, int starting_time_slot_no, int ending_time_slot_no, float t2, double vf, double kc, double s3_m, double vc,
           float est_speed[MAX_TIMEINTERVAL_PerDay], float est_volume_per_hour_per_lane[MAX_TIMEINTERVAL_PerDay])
    {

        // take nonnegative values
        volume = max(0.0, volume);

        if (VDF_type_no == 0) // BPR
        {
            return calculate_travel_time_based_on_BPR(volume);
        }

        if (VDF_type_no == 1) // queue-based VDF: QVDF, planning level, for the entire planning horizon. e.g. AM or PM
        {
            return calculate_travel_time_based_on_QVDF(volume, lanes, starting_time_slot_no, ending_time_slot_no, t2, vf, kc, s3_m, vc, est_speed, est_volume_per_hour_per_lane);
        }

        if (VDF_type_no == 2) // queue-based VDF: QVDF, planning level, for the entire planning horizon. e.g. AM or PM
        {
        }
        return FFTT;
    }

    double get_speed_from_volume(float hourly_volume, float vf, float kc, float s3_m)
    {
        //test data free_speed = 55.0f; 
        //speed = 52;
        //kc = 23.14167648;

        double max_lane_capacity = kc * vf / pow(2, 2 / s3_m);

        if (hourly_volume > max_lane_capacity)
            hourly_volume = max_lane_capacity;
        // we should add a capacity upper bound on hourly_volume;

        double coef_a = pow(kc, s3_m);
        double coef_b = pow(kc, s3_m) * pow(vf, s3_m / 2.0);
        double coef_c = pow(hourly_volume, s3_m);  // D is hourly demand volume, which is equivalent to flow q in S3 model

        double speed = pow((coef_b + pow(coef_b * coef_b - 4.0 * coef_a * coef_c, 0.5)) / (2.0 * coef_a), 2.0 / s3_m);    //under uncongested condition
        if (speed >= vf)
            speed = vf;

        if (speed < 0)
            speed = 0;

        return speed;

    }

    double get_volume_from_speed(float speed, float vf, float kc, float s3_m)
    {
        //test data free_speed = 55.0f; 
        //speed = 52;
        //kc = 23.14167648;

        if (speed < 0)
            return -1;

        double speed_ratio = vf / max(1, speed);
        if (speed_ratio <= 1.00001)
            speed_ratio = 1.00001;

        /*   float volume = 0;*/
        double ratio_difference = pow(speed_ratio, s3_m / 2) - 1;

        double ratio_difference_final = max(ratio_difference, 0.00000001);

        double volume = speed * kc * pow(ratio_difference_final, 1 / s3_m);  // volume per hour per lane

        return volume;

    }

    double calculate_travel_time_based_on_BPR(double volume)
    {

        // take nonnegative values
        volume = max(0.0, volume);
            VOC = volume / max(0.001, peak_load_factor) / max(0.00001, period_capacity);
            avg_travel_time = FFTT + FFTT * alpha * pow((volume + preload) / max(0.00001, period_capacity), beta);
            total_travel_time = (volume + preload) * avg_travel_time;
            marginal_base = FFTT * alpha * beta * pow((volume + preload) / max(0.00001, period_capacity), beta - 1);

            if (cycle_link_distance_in_km > 1)
            {
                if (volume > 10)
                    avg_travel_time += cycle_link_distance_in_km / 60.0 / 2.0; // add additional random delay due to signalized intersection by 
            }
            return avg_travel_time;
            // volume --> avg_traveltime
    }

    double calculate_travel_time_based_on_QVDF(double volume, int lanes, int starting_time_slot_no, int ending_time_slot_no, float t2, double vf, double kc, double s3_m, double vc,
        float est_speed[MAX_TIMEINTERVAL_PerDay], float est_volume_per_hour_per_lane[MAX_TIMEINTERVAL_PerDay])    
    {

                // take nonnegative values
            volume = max(0.0, volume);

            VOC = volume / max(0.001, peak_load_factor) / max(0.00001, period_capacity);

            float assign_period_start_time_in_hour = starting_time_slot_no* MIN_PER_TIMESLOT/60;
            float assign_period_end_time_in_hour = ending_time_slot_no * MIN_PER_TIMESLOT / 60;
            float L = assign_period_end_time_in_hour - assign_period_start_time_in_hour;
            
            float V = volume;

            float D = V / max(0.001, peak_load_factor); //under uncongested condition, D is the peak hour volume rate

            double coef_a = pow(kc, s3_m);
            double coef_b = pow(kc, s3_m) * pow(vf, s3_m / 2.0);
            double coef_c = pow(D, s3_m);  // D is hourly demand volume, which is equivalent to flow q in S3 model

            double avg_speed = pow((coef_b + pow(coef_b * coef_b - 4.0 * coef_a * coef_c, 0.5)) / (2.0 * coef_a), 2.0 / s3_m);    //under uncongested condition
            avg_travel_time = FFTT * (vf / max(0.0001, avg_speed));
           // work on congested condition
           double DCRatio = D/period_capacity;  // inflow demand to period_capacity ratio
           double P = pow(DCRatio, n);
           double mu = pow(DCRatio, n * (-1)) * D;
           double vt2 = vc / pow(P, s); //vt2=vc/P^s   

           double t0 = t2 - 0.5 * P;
           double t3 = t2 + 0.5 * P;

           double nonpeak_hourly_flow = 0;
               
               if(L - P >= 10.0 / 60.0)
               {
                   nonpeak_hourly_flow = (V - D) / max(1, lanes) / max(0.0001, L - P- 5.0/60.0);  //5.0/60.0 as one 5 min interval, as P includes both boundary points
               }

           dtalog.output() << "nonpeak_hourly_flow = " << nonpeak_hourly_flow << endl;

           // setup the upper bound on nonpeak flow rates
           double max_lane_capacity = kc * vf / pow(2, 2 / s3_m);

           if (nonpeak_hourly_flow > max_lane_capacity)
               nonpeak_hourly_flow = max_lane_capacity;

           double nonpeak_avg_speed = get_speed_from_volume(nonpeak_hourly_flow, vf, kc, s3_m);


           if (P < 1) // for P < 1, we will skip the creation of curve below vc. 
           {
               
               for (int t_in_min = assign_period_start_time_in_hour * 60; t_in_min < assign_period_end_time_in_hour * 60; t_in_min += 5)
               {
                   int t_interval = t_in_min / 5;

                   double t = t_in_min / 60.0;  // t in hour
                   double td_speed = 0;
                   double td_flow = nonpeak_hourly_flow;  // default

                   if (t0 <= t && t <= t3)
                   {
                       est_speed[t_interval] = avg_speed;
                       est_volume_per_hour_per_lane[t_interval] = get_volume_from_speed(avg_speed, vf, kc, s3_m);
                   }else if (t < t0)
                   {
                       td_speed = nonpeak_avg_speed;
                       est_volume_per_hour_per_lane[t_interval] = nonpeak_hourly_flow;
                   }
                   else if (t > t3)
                   {
                       td_speed = nonpeak_avg_speed;
                       est_volume_per_hour_per_lane[t_interval] = nonpeak_hourly_flow;
                   }

              }
               return avg_travel_time;
           }

           double q_t2 = 1.0 / vt2 * mu;
           double gamma = q_t2 * 4 / pow(P / 2, 4);  // because q_tw = w*mu =1/4 * gamma (P/2)^4, => 1/vt2 * mu = 1/4 * gamma  * (P/2)^4


           for (int t_in_min = assign_period_start_time_in_hour * 60; t_in_min < assign_period_end_time_in_hour * 60; t_in_min += 5)
           {
               double t = t_in_min / 60.0;  // t in hour
               double td_queue = 0;

               if (t0 <= t && t <= t3)
                 td_queue = 0.25 * gamma * pow((t - t0), 2) * pow((t - t3), 2);
           else
                 td_queue = 0;


                 double td_speed = 0;

                        if (P < 0.25)  // only for less than 1/4 hours
                        {
                            if (t < t0)
                            {
                                td_speed = nonpeak_avg_speed;
                            }
                            else if (t > t3)
                            {
                                td_speed = nonpeak_avg_speed;
                            }
                            else
                            {
                                td_speed = vc;
                            }
                        }
                        else // P > 0.25
                        {
                            if (t < t0)
                            {
                                td_speed = nonpeak_avg_speed;
                            }
                            else if (t > t3)
                            { 
                                td_speed = nonpeak_avg_speed;
                            }
                            else
                            {
                                td_speed = 1 / ( (td_queue / mu) + (1 / vc) );
                                
                            }

                            //                cout << "td_queue t" << t << " =  " << td_queue << ", speed =" << td_speed << endl;
                        }

                     int t_interval = t_in_min / 5;
                     double td_flow = get_volume_from_speed(td_speed, vf, kc, s3_m);
                     est_speed[t_interval] = td_speed;
                     est_volume_per_hour_per_lane[t_interval] = td_flow;
           }

           // peak load duration
           double pl_t0 = t2 - max(0.5, 0.5 * P);
           double pl_t3 = t2 + max(0.5, 0.5 * P);
           double est_peak_load_demand = 0;
           //est_non_peak_load_demand should not be counted, as avg non-peak rates have been determined by (V-D)/(L-P)

           double hourly_rate_2_volume_factor = lanes / 12.0;  // /12 to convert hourly to 5 min volume;
           // step 2
           for (int t_in_min = assign_period_start_time_in_hour * 60; t_in_min < assign_period_end_time_in_hour * 60; t_in_min += 5)
           {
               double t = t_in_min / 60.0;  // t in hour
               int t_interval = t_in_min / 5;
               
                if (t >= pl_t0  && t <= pl_t3)
                   {
                    est_peak_load_demand += est_volume_per_hour_per_lane[t_interval] * hourly_rate_2_volume_factor;
                   }
           }
           // step 3:
           double peak_load_volume_scale_factor = D / max(0.0001,est_peak_load_demand);


           //step 4
           for (int t_in_min = assign_period_start_time_in_hour * 60; t_in_min < assign_period_end_time_in_hour * 60; t_in_min += 5)
           {
               double t = t_in_min / 60.0;  // t in hour
               int t_interval = t_in_min / 5;

               if (t < pl_t0)
               {
                   est_volume_per_hour_per_lane[t_interval] = min(max_lane_capacity, est_volume_per_hour_per_lane[t_interval]);
               }
               else if (t > pl_t3)
               {
                   est_volume_per_hour_per_lane[t_interval] = min(max_lane_capacity, est_volume_per_hour_per_lane[t_interval]);
               }
               else
               {
                   est_volume_per_hour_per_lane[t_interval] = min(max_lane_capacity, est_volume_per_hour_per_lane[t_interval] * peak_load_volume_scale_factor);
               }
           }


           //final stage: avg delay in peak load hours

                float alpha = 8.0 / 15.0;
                float beta = n * s;
                avg_speed = vc / (1 + alpha * (pow(DCRatio, beta) - 1));
                
                avg_travel_time = FFTT * (vf / max(0.0001, avg_speed));

            return avg_travel_time;
     }


    double m;
    // we should also pass uncongested_travel_time as link_distance_in_km/(speed_at_capacity)
    double VOC;
    //updated BPR-X parameters
    double gamma;
    double mu;
    //peak hour factor
    double PHF;
    //standard BPR parameter
    double alpha;
    double beta;

    double peak_load_factor;  // peak load factor
    double sav; // volume used in sensitivty analysis
    //https://en.wikipedia.org/wiki/Peak_demand
    //https://aquicore.com/blog/what-is-peak-load/ 
    //peak load pricing: https://saylordotorg.github.io/text_introduction-to-economic-analysis/s16-07-peak-load-pricing.html
    //peak hour factor: https://help.miovision.com/s/article/Peak-Hour-Factor-PHF-Explained?language=en_US 
    /*The Peak Hour Factor(PHF) compares the traffic volume during the busiest 15 - minutes of the peak hour with the total volume during the peak hour.It indicates how consistent traffic volume is during the peak hour.*/

    double n;
    double s;
    double vc;

    double preload;
    double toll[MAX_AGNETTYPES];
    double pce[MAX_AGNETTYPES];
    double penalty;
    double LR_price[MAX_AGNETTYPES];
    double LR_RT_price[MAX_AGNETTYPES];;
    string allowed_uses;

    double rho;
    double marginal_base;
    // in 15 min slot
    int starting_time_slot_no;
    int ending_time_slot_no;
    float cycle_link_distance_in_km;
    float red_time;
    float effective_green_time;
    int start_green_time;
    int end_green_time;

    int t0, t3;

    bool bValidQueueData;
    string period;

    double period_capacity;  // link based period_capacity
    double time_period_in_hour;
    double FFTT;

    double congestion_period_P;
    // inpput
    double volume;

    //output
    double avg_delay;
    double avg_travel_time;
    double total_travel_time;

    //// t starting from starting_time_slot_no if we map back to the 24 hour horizon
    float queue_length;
    float arrival_flow_volume;
    float discharge_rate;  // period based
    float avg_waiting_time;
    float travel_time;

    std::map<int, float> travel_time_per_iteration_map;
};



