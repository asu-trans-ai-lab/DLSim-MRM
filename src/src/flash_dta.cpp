/* Portions Copyright 2019 Xuesong Zhou and Peiheng Li
 *
 * If you help write or modify the code, please also list your names here.
 * The reason of having Copyright info here is to ensure all the modified version, as a whole, under the GPL
 * and further prevent a violation of the GPL.
 *
 * More about "How to use GNU licenses for your own software"
 * http://www.gnu.org/licenses/gpl-howto.html
 */

#ifdef _WIN32
#include "pch.h"
#endif

#include <iostream>
#include <fstream>
#include <iomanip>
#include "config.h"
#include "utils.h"
using namespace std;
void write_default_setting_file_if_not_exist()
{

    CCSVParser parser_settings;
    parser_settings.IsFirstLineHeader = false;

    if (parser_settings.OpenCSVFile("settings.csv", false))
    {
        parser_settings.CloseCSVFile();
        return;
    }

  ofstream myfile;
  myfile.open("settings.csv");
  myfile << "[assignment] ,,assignment_mode,column_generation_iterations,column_updating_iterations,odme_iterations,sensitivity_analysis_iterations,simulation_iterations,remarks" << endl;
  myfile << ",,dta,5,0,0,0,1," << endl;
  myfile << "[agent_type],agent_type,name,,vot,flow_type,pce,occ," << endl;
  myfile << ",auto,auto,,10,0,1,1," << endl;
  myfile << ",walk,walk,,10,0,1,1," << endl;
  myfile << ",bike,bike,,10,0,1,1," << endl;
  myfile << "" << endl;
  myfile << "[link_type],link_type,link_type_name,,agent_type_blocklist,type_code,traffic_flow_code,vdf_type," << endl;
  myfile << ",1,Highway / Expressway,,,f,0,bpr," << endl;
  myfile << ",2,Major arterial,,,a,0,bpr," << endl;
  myfile << "[demand_period], demand_period_id, demand_period, , time_period, peak_time, , , time_interval0, time_interval2, time_interval3, time_interval4, time_interval5, time_interval6, time_interval7, time_interval8, time_interval9, time_interval10, time_interval11, time_interval12, time_interval13, time_interval14, time_interval15, time_interval16, time_interval17, time_interval18, time_interval19, time_interval20, time_interval21, time_interval22, time_interval23, time_interval24, time_interval25, time_interval26, time_interval27, time_interval28, time_interval29, time_interval30, time_interval31, time_interval32, time_interval33, time_interval34, time_interval35, time_interval36, time_interval37, time_interval38, time_interval39, time_interval40, time_interval41, time_interval42, time_interval43, time_interval44, time_interval45, time_interval46, time_interval47, time_interval48, time_interval49, time_interval50, time_interval51, time_interval52, time_interval53, time_interval54, time_interval55, time_interval56, time_interval57, time_interval58, time_interval59, time_interval60, time_interval61, time_interval62, time_interval63, time_interval64, time_interval65, time_interval66, time_interval67, time_interval68, time_interval69, time_interval70, time_interval71, time_interval72, time_interval73, time_interval74, time_interval75, time_interval76, time_interval77, time_interval78, time_interval79, time_interval80, time_interval81, time_interval82, time_interval83, time_interval84, time_interval85, time_interval86, time_interval87, time_interval88, time_interval89, time_interval90, time_interval91, time_interval92, time_interval93, time_interval94, time_interval95, time_interval96" << endl;
  myfile << ", 1,AM, ,0700_1000,, , , 0.001712541, 0.001518821, 0.001335712, 0.001173243, 0.001071592, 0.000984901, 0.000957369, 0.000905166, 0.000875259, 0.000889426, 0.000929134, 0.000931492, 0.000955705, 0.001149788, 0.001487502, 0.001704861, 0.001967319, 0.0028493, 0.00410496, 0.004761717, 0.005249334, 0.006865023, 0.008763547, 0.009727353, 0.0096528, 0.011408503, 0.013377631, 0.015006034, 0.015620492, 0.017031255, 0.017982062, 0.01805379, 0.016524355, 0.015870414, 0.015175223, 0.014499898, 0.013261577, 0.01297958, 0.013092157, 0.013029354, 0.01241647, 0.012602333, 0.012872216, 0.013050637, 0.013227912, 0.013698897, 0.014021648, 0.014281785, 0.014480951, 0.014645749, 0.014698547, 0.014661224, 0.014505402, 0.014697936, 0.015068191, 0.015193797, 0.015487111, 0.016308453, 0.017314827, 0.017722353, 0.017630525, 0.018148853, 0.018588326, 0.018744003, 0.01892459, 0.019213207, 0.019172793, 0.019202233, 0.01957678, 0.019723407, 0.018811942, 0.017812392, 0.016733696, 0.015880365, 0.014502684, 0.013161254, 0.012091192, 0.011245231, 0.010145693, 0.009361706, 0.008888914, 0.008671416, 0.008012391, 0.007403109, 0.007094947, 0.006746031, 0.006045864, 0.005352276, 0.004919291, 0.004421213, 0.003936383, 0.003397222, 0.003016297, 0.002667224, 0.002334765, 0.002029159" << endl;
  myfile << "[demand_file_list],file_sequence_no,file_name,,format_type,demand_period,agent_type,scale_factor," << endl;
  myfile << ",1,input_matrix.csv,,matrix,AM,auto,0.333," << endl;
  myfile << ",2,input_matrix.csv,,matrix,AM,walk,0.333," << endl;
  myfile << ",3,input_matrix.csv,,matrix,AM,bike,0.333," << endl;

  myfile.close();
  return;
}

int main()
{
    // reset all the log files to defult 0: not output; if want to output these logs set to 1
    dtalog.output() << "DTALite Log" << std::fixed << std::setw(12) << '\n';
    dtalog.debug_level() = 0;
    dtalog.log_sig() = 0;
    dtalog.log_odme() = 0;
    dtalog.log_path() = 0;
    dtalog.log_dta() = 0;
    dtalog.log_ue() = 0;

    int column_generation_iterations = 20;
    int column_updating_iterations = 40;
    int ODME_iterations = 20;
    int sensitivity_analysis_iterations = 0;
	int number_of_memory_blocks = 8;
    int simulation_iterations = 0;

    int signal_updating_output = 0;
    // generate link performance and agent file
    int assignment_mode = 1;
    bool flag_default = false;
    int default_volume = 1;
    int link_length_in_meter_flag = 0;

    write_default_setting_file_if_not_exist();
    CCSVParser parser_settings;
    parser_settings.IsFirstLineHeader = false;

    if (parser_settings.OpenCSVFile("settings.csv", false))
    {
        while (parser_settings.ReadRecord_Section())
        {
            if (parser_settings.SectionName.compare("assignment") >= 0);
            {
                std::string assignment_mode_str;
                parser_settings.GetValueByFieldName("column_generation_iterations", column_generation_iterations, true, true);
                parser_settings.GetValueByFieldName("assignment_mode", assignment_mode_str);
                parser_settings.GetValueByFieldName("assignment_mode", link_length_in_meter_flag);
                
                // these are the assignment modes
                // two usually methods are ue (user equilibrium) and dta (dynamic traffic assignment)
                // the main difference of these two methods are different output in link_performance.csv
                // for basic uses set assignment mode to 'ue'
                // for more detailed link performances (one minute) set 'dta'
                if (assignment_mode_str == "lue")
                    assignment_mode = 0;
                else if (assignment_mode_str == "dta")
                    assignment_mode = 3;
                else if (assignment_mode_str == "zone2access")  //access link
                    assignment_mode = 21;
                else if (assignment_mode_str == "cbi")  //congestion bottleneck identification
                    assignment_mode = 11;
                else if (assignment_mode_str == "cbsa")  //congestion bottleneck sensitivity
                    assignment_mode = 12;
                else
                {
                    dtalog.output() << "assignment_mode " << assignment_mode_str.c_str() << " in settings.csv is invalid." << std::endl;
                    g_program_stop();
                }

                // iteration number of reassignment
                parser_settings.GetValueByFieldName("column_updating_iterations", column_updating_iterations, true, true);

                if(assignment_mode == 3)
                {
                    parser_settings.GetValueByFieldName("odme_iterations", ODME_iterations, true, false);
                    parser_settings.GetValueByFieldName("sensitivity_analysis_iterations", sensitivity_analysis_iterations, true, false);
                }
                parser_settings.GetValueByFieldName("simulation_iterations", simulation_iterations, true, false);
                
                    // the start interation of generating signals, if there is no signals set this number larger than the iteration number
                parser_settings.GetValueByFieldName("number_of_memory_blocks", number_of_memory_blocks, false, false);
				dtalog.output() << "number_of_memory_blocks = " << number_of_memory_blocks << " in settings.csv." << std::endl;


                // just one record
                break;
            }

            if (parser_settings.SectionName == "[log]")
            {
                parser_settings.GetValueByFieldName("sig", dtalog.log_sig(), false);
                parser_settings.GetValueByFieldName("odme", dtalog.log_odme(), false);
                parser_settings.GetValueByFieldName("path", dtalog.log_path(), false);
                parser_settings.GetValueByFieldName("ue", dtalog.log_ue(), false);

                // just one record
                break;
            }
        }
    }
    // obtain initial flow values
    network_assignment(assignment_mode, column_generation_iterations, column_updating_iterations, ODME_iterations, sensitivity_analysis_iterations, simulation_iterations, number_of_memory_blocks);

    return 0;
}