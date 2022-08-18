////  Portions Copyright 2019
////   Xuesong (Simon) Zhou
////   If you help write or modify the code, please also list your names here.
////   The reason of having Copyright info here is to ensure all the modified version, as a whole, under the GPL 
////   and further prevent a violation of the GPL.
//
//// More about "How to use GNU licenses for your own software"
//// http://www.gnu.org/licenses/gpl-howto.html
//
//#pragma warning( disable : 4305 4267 4018) 
//
//#include <iostream>
//#include <fstream>
//#include <list> 
//#include <omp.h>
//#include <algorithm>
//#include <time.h>
//#include <functional>
//#include <stdio.h>   
//#include <math.h>
//
//
//#include <stack>
//#include <string>
//#include <vector>
//#include <map>
//#include <sstream>
//#include <iostream>
//#include <iomanip>
//using namespace std;
//using std::string;
//using std::ifstream;
//using std::vector;
//using std::map;
//using std::istringstream;
//using std::max;
//
//string g_info_String = "";
//class CCSVParser_signal
//{
//public:
//	char Delimiter;
//	bool IsFirstLineHeader;
//	ifstream inFile;
//	string mFileName;
//	vector<string> LineFieldsValue;
//	vector<string> Headers;
//	map<string, int> FieldsIndices;
//
//	vector<int> LineIntegerVector;
//
//public:
//	void  ConvertLineStringValueToIntegers()
//	{
//		LineIntegerVector.clear();
//		for (unsigned i = 0; i < LineFieldsValue.size(); i++)
//		{
//			std::string si = LineFieldsValue[i];
//			int value = atoi(si.c_str());
//
//			if (value >= 1)
//				LineIntegerVector.push_back(value);
//
//		}
//	}
//	vector<string> GetHeaderVector()
//	{
//		return Headers;
//	}
//
//	bool m_bDataHubSingleCSVFile;
//	string m_DataHubSectionName;
//	bool m_bLastSectionRead;
//
//	bool m_bSkipFirstLine;  // for DataHub CSV files
//
//	CCSVParser_signal(void)
//	{
//		Delimiter = ',';
//		IsFirstLineHeader = true;
//		m_bSkipFirstLine = false;
//		m_bDataHubSingleCSVFile = false;
//		m_bLastSectionRead = false;
//	}
//
//	~CCSVParser_signal(void)
//	{
//		if (inFile.is_open()) inFile.close();
//	}
//
//
//	bool OpenCSVFile(string fileName, bool b_required)
//	{
//		mFileName = fileName;
//		inFile.open(fileName.c_str());
//
//		if (inFile.is_open())
//		{
//			if (IsFirstLineHeader)
//			{
//				string s;
//				std::getline(inFile, s);
//				vector<string> FieldNames = ParseLine(s);
//
//				for (size_t i = 0; i < FieldNames.size(); i++)
//				{
//					string tmp_str = FieldNames.at(i);
//					size_t start = tmp_str.find_first_not_of(" ");
//
//					string name;
//					if (start == string::npos)
//					{
//						name = "";
//					}
//					else
//					{
//						name = tmp_str.substr(start);
//						//			TRACE("%s,", name.c_str());
//					}
//
//
//					FieldsIndices[name] = (int)i;
//				}
//			}
//
//			return true;
//		}
//		else
//		{
//			if (b_required)
//			{
//
//				cout << "File " << fileName << " does not exist. Please check." << endl;
//				//g_ProgramStop();
//			}
//			return false;
//		}
//	}
//
//
//	void CloseCSVFile(void)
//	{
//		inFile.close();
//	}
//
//
//
//	bool ReadRecord()
//	{
//		LineFieldsValue.clear();
//
//		if (inFile.is_open())
//		{
//			string s;
//			std::getline(inFile, s);
//			if (s.length() > 0)
//			{
//
//				LineFieldsValue = ParseLine(s);
//
//				return true;
//			}
//			else
//			{
//
//				return false;
//			}
//		}
//		else
//		{
//			return false;
//		}
//	}
//
//	vector<string> ParseLine(string line)
//	{
//		vector<string> SeperatedStrings;
//		string subStr;
//
//		if (line.length() == 0)
//			return SeperatedStrings;
//
//		istringstream ss(line);
//
//
//		if (line.find_first_of('"') == string::npos)
//		{
//
//			while (std::getline(ss, subStr, Delimiter))
//			{
//				SeperatedStrings.push_back(subStr);
//			}
//
//			if (line.at(line.length() - 1) == ',')
//			{
//				SeperatedStrings.push_back("");
//			}
//		}
//		else
//		{
//			while (line.length() > 0)
//			{
//				size_t n1 = line.find_first_of(',');
//				size_t n2 = line.find_first_of('"');
//
//				if (n1 == string::npos && n2 == string::npos) //last field without double quotes
//				{
//					subStr = line;
//					SeperatedStrings.push_back(subStr);
//					break;
//				}
//
//				if (n1 == string::npos && n2 != string::npos) //last field with double quotes
//				{
//					size_t n3 = line.find_first_of('"', n2 + 1); // second double quote
//
//					//extract content from double quotes
//					subStr = line.substr(n2 + 1, n3 - n2 - 1);
//					SeperatedStrings.push_back(subStr);
//
//					break;
//				}
//
//				if (n1 != string::npos && (n1 < n2 || n2 == string::npos))
//				{
//					subStr = line.substr(0, n1);
//					SeperatedStrings.push_back(subStr);
//					if (n1 < line.length() - 1)
//					{
//						line = line.substr(n1 + 1);
//					}
//					else // comma is the last char in the line string, push an empty string to the back of vector
//					{
//						SeperatedStrings.push_back("");
//						break;
//					}
//				}
//
//				if (n1 != string::npos && n2 != string::npos && n2 < n1)
//				{
//					size_t n3 = line.find_first_of('"', n2 + 1); // second double quote
//					subStr = line.substr(n2 + 1, n3 - n2 - 1);
//					SeperatedStrings.push_back(subStr);
//					size_t idx = line.find_first_of(',', n3 + 1);
//
//					if (idx != string::npos)
//					{
//						line = line.substr(idx + 1);
//					}
//					else
//					{
//						break;
//					}
//				}
//			}
//
//		}
//
//		return SeperatedStrings;
//	}
//
//	template <class T> bool GetValueByFieldName(string field_name, T& value, bool NonnegativeFlag = true, bool required_field = true)
//	{
//
//		if (FieldsIndices.find(field_name) == FieldsIndices.end())
//		{
//			if (required_field)
//			{
//				cout << "Field " << field_name << " in file " << mFileName << " does not exist. Please check the file." << endl;
//
//				exit(0);
//			}
//			return false;
//		}
//		else
//		{
//			if (LineFieldsValue.size() == 0)
//			{
//				return false;
//			}
//
//			int size = (int)(LineFieldsValue.size());
//			if (FieldsIndices[field_name] >= size)
//			{
//				return false;
//			}
//
//			string str_value = LineFieldsValue[FieldsIndices[field_name]];
//
//			if (str_value.length() <= 0)
//			{
//				return false;
//			}
//
//			istringstream ss(str_value);
//
//			T converted_value;
//			ss >> converted_value;
//
//			if (/*!ss.eof() || */ ss.fail())
//			{
//				return false;
//			}
//
//			if (NonnegativeFlag && converted_value < 0)
//				converted_value = 0;
//
//			value = converted_value;
//			return true;
//		}
//	}
//
//
//	bool GetValueByFieldName(string field_name, string& value)
//	{
//		if (FieldsIndices.find(field_name) == FieldsIndices.end())
//		{
//			return false;
//		}
//		else
//		{
//			if (LineFieldsValue.size() == 0)
//			{
//				return false;
//			}
//
//			unsigned int index = FieldsIndices[field_name];
//			if (index >= LineFieldsValue.size())
//			{
//				return false;
//			}
//			string str_value = LineFieldsValue[index];
//
//			if (str_value.length() <= 0)
//			{
//				return false;
//			}
//
//			value = str_value;
//			return true;
//		}
//	}
//
//};
//class C_SignalMainModual {
//public:
//	C_SignalMainModual()
//	{
//
//		g_number_of_links = 0;
//		g_number_of_service_arcs = 0;
//		g_number_of_nodes = 0;
//
//		b_debug_detail_flag = 1;
//
//		g_pFileDebugLog = NULL;
//		g_informationCount = 0;
//
//
//		SYSTEMTIME st;
//		GetLocalTime(&st);
//		g_pFileDebugLog = fopen("log.txt", "at");
//		fprintf(g_pFileDebugLog, "Logging Time: %d/%d/%d_%d:%d:%d:%d---------------------------------------¡ý\n",st.wYear,st.wMonth,st.wDay, st.wHour, st.wMinute, st.wSecond, st.wMilliseconds);
//		fclose(g_pFileDebugLog);
//		messageType_int_to_string[0] = "MSG";//message
//		messageType_int_to_string[1] = "OUT";//output
//		messageType_int_to_string[2] = "WAR";//Warning
//
//	}
//
//	int b_debug_detail_flag;
//	std::map<int, int> g_internal_node_to_seq_no_map;  // hash table, map external node number to internal node sequence no. 
//
//
//	std::map<string, int> g_road_link_id_map;
//
//
//	int g_number_of_links;
//	int g_number_of_service_arcs;
//	int g_number_of_nodes;
//	int g_number_of_zones;
//
//	int g_LoadingStartTimeInMin;
//	int g_LoadingEndTimeInMin;
//	int g_informationCount;
//
//	FILE* g_pFileDebugLog;
//	std::map<int, char*> messageType_int_to_string;
//
//	void WriteLog(int messageType, string info, int level)
//	{
//		string padding = "";
//		if (level > 1)
//		{
//			for (size_t l = 2; l <= level; l++)
//			{
//				padding.append("--");
//			}
//		}
//		if (level == -1)
//		{
//			padding = "++++++++++++++++++++++++++++++";
//		}
//		if (level == -2)
//		{
//			padding = "++++++++++++++++++++++++++++++";
//		}
//
//		g_pFileDebugLog = fopen("log.txt", "at");
//		fprintf(g_pFileDebugLog, "Info_ID: %d| Info_Type: %s \t| Info: %s%s\n",  g_informationCount, messageType_int_to_string[messageType], padding.c_str(), info.c_str());
//		fclose(g_pFileDebugLog);
//		cout << "Info_ID: " << g_informationCount << "| Info_Type: " << messageType_int_to_string[messageType] << " \t| Info: "<< padding.c_str()  << info.c_str() << endl;
//		//OutputDebugStringA(szLog);
//		g_informationCount++;
//	}
//
//	string Output_Combine_Message_And_Value_Int(string info,int value)
//	{
//		string value_String = to_string(value);
//
//		info.append(": ");
//
//		info.append(value_String);
//
//		return info;
//	}
//	string Output_Combine_Message_And_Value_Double(string info, double value)
//	{
//		string value_String = to_string(value);
//
//		info.append(": ");
//
//		info.append(value_String);
//
//		return info;
//	}
//
//};
//class C_SignalLink
//{
//public:
//	int zone_seq_no_for_outgoing_connector;
//	int m_LeftTurn_link_seq_no;
//	float greenTime;
//
//	int m_RandomSeed;
//	int link_seq_no;
//	string link_id;
//	int traffic_flow_code;
//	int spatial_capacity_in_vehicles;
//	int from_node_seq_no;
//	int to_node_seq_no;
//	int link_type;
//	bool service_arc_flag;
//	float toll;
//	float route_choice_cost;
//
//
//	float free_flow_travel_time_in_min;
//	float lane_capacity;
//	int number_of_lanes;
//	float length;
//public:
//	C_SignalLink()  // construction 
//	{
//
//	}
//
//	~C_SignalLink()
//	{
//		//if (flow_volume_for_each_o != NULL)
//		//	delete flow_volume_for_each_o;
//	}
//
//
//};
//
//C_SignalMainModual MainModual;
//struct SMovementData
//{
//	//SMovementData()
//	//{
//
//	//}
//	bool Enable = false;//no value->false
//	int LinkSeqNo;
//	int PhaseNo;
//	//enum stage StageNo;
//	vector<enum stage> StageNo_in_Order; 
//	int Volume_per_lane;
//	int Lanes;
//	int SharedLanes;
//	enum group GroupNo;//L->1;T/R->2
//	enum direction DirectionNo;
//
//	double signal_link_capacity; double saturation_flow_rate; int green_time; int pre_set_stage_no=-1; int cycle_length; /*int priority;*/
//
//	//by Acquisition through processing
//	int Assignment_Order;// not in use, Differentiate primary and secondary movements, then determine stages  according to volumes
//	enum left_Turn_Treatment Left_Turn_Treatment;
//
//	string linkID;
//
//};
//
//#pragma region Fields
//
////enum
//enum movement_Index { EBL = 1, EBT = 2, EBR = 3, WBL = 4, WBT = 5, WBR = 6, NBL = 7, NBT = 8, NBR = 9, SBL = 10, SBT = 11, SBR = 12 };
//enum stage { no_stage = -1, stage1 = 1, stage2 = 2, stage3 = 3, stage4 = 4 };
//enum direction { E = 1, W = 2, N = 3, S = 4 };
//enum group { L = 1, T_AND_R = 2 };
//enum left_Turn_Treatment { perm = 0, prot = 1, no_business = -1 };
//
////array size, for constructing matirx or array
//const int laneColumnSize = 32;
//const int movementSize = 32;
//const int NEMA_PhaseSize = 32;//temp enable=false
//const int stageSize = 5;
//const int ringSize = 3;
//const int directionSize = 5;
//const int groupSize = 5;
//
////index range, for indexing, from 1 to 0+range (including the 0+range th)
//
////parameters
//double l = 12;
//double x_c_Input = 0.9;
//double PHF = 1;
//
//double f_1 = 1;
//double f_2 = 1;
//double t_L = 4;
//double t_Yellow = 4;
//double t_AR = 2;
//
//double minGreenTime = 5;
//#pragma endregion
//
//class C_SignalNode
//{
//public:
//	int node_seq_no;  // sequence number 
//	int node_id;      //external node number 
//	int zone_id = -1;
//
//	double x;
//	double y;
//
//	std::vector<int> m_outgoing_link_seq_no_vector;
//	std::vector<int> m_incoming_link_seq_no_vector;
//
//	std::vector<int> m_to_node_seq_no_vector;
//	std::map<int, int> m_to_node_2_link_seq_no_map;
//
//public:
//	C_SignalNode()//Constructor
//	{
//
//		zone_id = -1;
//		node_seq_no = -1;
//		y_StageMax = 1;
//		x_c_output = 0.9;
//		c_Min = 60;
//
//		movement_str_to_index_map["EBL"] = EBL;
//		movement_str_to_index_map["EBT"] = EBT;
//		movement_str_to_index_map["EBR"] = EBR;
//
//		movement_str_to_index_map["WBL"] = WBL;
//		movement_str_to_index_map["WBT"] = WBT;
//		movement_str_to_index_map["WBR"] = WBR;
//
//		movement_str_to_index_map["NBL"] = NBL;
//		movement_str_to_index_map["NBT"] = NBT;
//		movement_str_to_index_map["NBR"] = NBR;
//
//		movement_str_to_index_map["SBL"] = SBL;
//		movement_str_to_index_map["SBT"] = SBT;
//		movement_str_to_index_map["SBR"] = SBR;
//
//		movement_str_array[EBL] = "EBL";
//		movement_str_array[EBT] = "EBT";
//		movement_str_array[EBR] = "EBR";
//		movement_str_array[WBL] = "WBL";
//		movement_str_array[WBT] = "WBT";
//		movement_str_array[WBR] = "WBR";
//		movement_str_array[NBL] = "NBL";
//		movement_str_array[NBT] = "NBT";
//		movement_str_array[NBR] = "NBR";
//		movement_str_array[SBL] = "SBL";
//		movement_str_array[SBT] = "SBT";
//		movement_str_array[SBR] = "SBR";
//
//
//		movement_str_to_direction_map["EBL"] = E;
//		movement_str_to_direction_map["EBT"] = E;
//		movement_str_to_direction_map["EBR"] = E;
//
//		movement_str_to_direction_map["WBL"] = W;
//		movement_str_to_direction_map["WBT"] = W;
//		movement_str_to_direction_map["WBR"] = W;
//
//		movement_str_to_direction_map["NBL"] = N;
//		movement_str_to_direction_map["NBT"] = N;
//		movement_str_to_direction_map["NBR"] = N;
//
//		movement_str_to_direction_map["SBL"] = S;
//		movement_str_to_direction_map["SBT"] = S;
//		movement_str_to_direction_map["SBR"] = S;
//
//
//		left_Movement_Opposing_Index_Map[EBL] = WBT;
//		left_Movement_Opposing_Index_Map[WBL] = EBT;
//		left_Movement_Opposing_Index_Map[NBL] = SBT;
//		left_Movement_Opposing_Index_Map[SBL] = NBT;
//
//
//		left_Movement_Counterpart_Index_Map[EBL] = EBT;
//		left_Movement_Counterpart_Index_Map[WBL] = WBT;
//		left_Movement_Counterpart_Index_Map[NBL] = NBT;
//		left_Movement_Counterpart_Index_Map[SBL] = SBT;
//
//		left_Movement_Counterpart_Right_Trun_Index_Map[EBL] = EBR;
//		left_Movement_Counterpart_Right_Trun_Index_Map[WBL] = WBR;
//		left_Movement_Counterpart_Right_Trun_Index_Map[NBL] = NBR;
//		left_Movement_Counterpart_Right_Trun_Index_Map[SBL] = SBR;
//
//		movement_Index_to_Group_Map[EBL] = group(1);
//		movement_Index_to_Group_Map[EBT] = group(2);
//		movement_Index_to_Group_Map[EBR] = group(2);
//		movement_Index_to_Group_Map[WBL] = group(1);
//		movement_Index_to_Group_Map[WBT] = group(2);
//		movement_Index_to_Group_Map[WBR] = group(2);
//		movement_Index_to_Group_Map[NBL] = group(1);
//		movement_Index_to_Group_Map[NBT] = group(2);
//		movement_Index_to_Group_Map[NBR] = group(2);
//		movement_Index_to_Group_Map[SBL] = group(1);
//		movement_Index_to_Group_Map[SBT] = group(2);
//		movement_Index_to_Group_Map[SBR] = group(2);
//
//		direction_index_to_str_map[E] = "E";
//		direction_index_to_str_map[W] = "W";
//		direction_index_to_str_map[N] = "N";
//		direction_index_to_str_map[S] = "S";
//
//		intersection_Average_Delay = 0;
//		intersection_Total_Delay = 0;
//		intersection_Total_Volume = 0;
//
//		left_Turn_Treatment_index_to_str_map[prot] = "Protected";
//		left_Turn_Treatment_index_to_str_map[perm] = "Permissive";
//		left_Turn_Treatment_index_to_str_map[no_business] = "Null";
//
//
//		for (int s = 0; s <= stageSize; s++)
//		{
//			y_Max_Stage_Array[s] = 0;
//
//			for (int m = 0; m <= movementSize; m++)
//			{
//				saturation_Flow_Rate_Matrix[s][m] = 0;
//				capacity_by_Stage_and_Movement_Matrix[s][m] = 0;
//			}
//
//			for (int d = 0; d <= directionSize; d++)
//			{
//				for (int g = 0; g <= groupSize; g++)
//				{
//					stage_Direction_Candidates_Matrix[s][d][g] = 0;
//
//				}
//				approach_Average_Delay_Array[d] = 0;
//				approach_Total_Delay_Array[d] = 0;
//				approach_Total_Volume_Array[d] = 0;
//			}
//
//		}
//
//		for (int m = 0; m <= movementSize; m++)
//		{
//			movement_Array[m].Enable = false;
//			movement_Array[m].Volume_per_lane = 0;
//		}
//	}
//
//	int movement_Range = 12;
//	int direction_Range = 4;
//
//	std::map<string, enum movement_Index> movement_str_to_index_map;
//	std::map<enum direction, string > direction_index_to_str_map;
//	std::map<enum left_Turn_Treatment, string > left_Turn_Treatment_index_to_str_map;
//
//
//	string movement_str_array[movementSize + 1];
//		
//	std::map<string, enum direction> movement_str_to_direction_map;
//
//	std::map<enum movement_Index, enum movement_Index> left_Movement_Opposing_Index_Map;
//
//	std::map<enum movement_Index, enum movement_Index> left_Movement_Counterpart_Index_Map;//L -> T
//
//	std::map<enum movement_Index, enum movement_Index> left_Movement_Counterpart_Right_Trun_Index_Map;//L -> R
//
//
//	std::map<enum movement_Index, enum group> movement_Index_to_Group_Map;
//
//
//	//array
//
//	SMovementData movement_Array[movementSize + 1]; //indexed by enum movement_Direction   --Step 2
//	int green_Start_Stage_Array[stageSize + 1];
//	int green_End_Stage_Array[stageSize+1];
//	double y_Max_Stage_Array[stageSize+1];
//	double green_Time_Stage_Array[stageSize + 1];
//	double cumulative_Green_Start_Time_Stage_Array[stageSize + 1];
//	double cumulative_Green_End_Time_Stage_Array[stageSize + 1];
//	double effective_Green_Time_Stage_Array[stageSize + 1];
//	double cumulative_Effective_Green_Start_Time_Stage_Array[stageSize + 1];
//	double cumulative_Effective_Green_End_Time_Stage_Array[stageSize + 1];
//	double ratio_of_Effective_Green_Time_to_Cycle_Length_Array[stageSize + 1];
//
//	double approach_Average_Delay_Array[directionSize+1];
//	double approach_Total_Delay_Array[directionSize + 1];
//	double approach_Total_Volume_Array[directionSize + 1];
//
//
//
//	//matrix
//	double saturation_Flow_Rate_Matrix[stageSize+1][movementSize + 1];//s
//	double y_Stage_Movement_Matrix[stageSize + 1][movementSize + 1];
//	double stage_Direction_Candidates_Matrix[stageSize + 1][directionSize + 1][groupSize + 1];
//
//	double capacity_by_Stage_and_Movement_Matrix[stageSize + 1][movementSize + 1];
//	double v_over_C_by_Stage_and_Movement_Matrix[stageSize + 1][movementSize + 1];
//
//	double average_Uniform_Delay_Matrix[stageSize + 1][movementSize + 1];
//
//
//	int NEMA_Phase_Matrix[5][3] = { 0,0,0,0,1,5,0,2,6,0,3,7,0,4,8 };//row->NEMA_phases; col->rings
//	int green_Start_NEMA_Phase[5][3];
//
//	SMovementData Stage_Ring_Movement_Matrix[ringSize + 1][stageSize + 1];
//
//	//variables
//	int stage_Range;
//	double y_StageMax;
//	double x_c_output;
//	double c_Min;
//	double c_Optimal;
//
//	double intersection_Average_Delay;
//	double intersection_Total_Delay;
//	double intersection_Total_Volume;
//
//	string LOS;
//
//	void PerformQEM(int nodeID)
//	{
//		MainModual.WriteLog(0, "Data Loading : Movement Volume of Main Node", 1);
//		//AddMovementVolume_from_LinkPerformance();
//
//		MainModual.WriteLog(0, "Step 2: Perform QEM", 1);
//
//		MainModual.WriteLog(0, "Step 2.1: Set Left Turn Treatments",2);
//		Set_Left_Turn_Treatment();
//
//		MainModual.WriteLog(0, "Step 2.2: Set StageNos",2);
//		Set_StageNo_for_Movements(); 
//
//		MainModual.WriteLog(0, "Step 2.3: Set Saturation Flow Rate Matrix",2);
//		Set_Saturation_Flow_Rate_Matrix(); 
//
//		MainModual.WriteLog(0, "Step 2.4: Calculate Flow Ratio",2);
//		Calculate_Flow_of_Ratio_Max();
//		MainModual.WriteLog(1, MainModual.Output_Combine_Message_And_Value_Double("Max y Value", y_StageMax),3);
//
//		MainModual.WriteLog(0, "Step 2.5: Calculate Total Cycle Lost Time", 2);
//		Calculate_Total_Cycle_Lost_Time();
//		MainModual.WriteLog(1, MainModual.Output_Combine_Message_And_Value_Double("Total Cycle Lost Time", l), 3);
//
//		MainModual.WriteLog(0, "Step 2.6: Calculate the Minimum and Optimal Cycle Length",2);
//		Calculate_the_Minimum_And_Optimal_Cycle_Length();
//		MainModual.WriteLog(1, MainModual.Output_Combine_Message_And_Value_Double("Min Cycle Length", c_Min),3);
//		MainModual.WriteLog(1, MainModual.Output_Combine_Message_And_Value_Double("Optimal Cycle Length", c_Optimal), 3);
//
//
//		MainModual.WriteLog(0, "Step 2.7: Recalculate xc",2);
//		Calculate_the_x_c_Output();
//		MainModual.WriteLog(1, MainModual.Output_Combine_Message_And_Value_Double("Recalculated xc", x_c_output),3);
//
//		MainModual.WriteLog(0, "Step 2.8: Timing for Stages",2);
//		Calculate_Green_Time_for_Stages();
//		Printing_Green_Time_for_Stages();
//
//		MainModual.WriteLog(0, "Step 2.9: Calculate capacity and ratio V/C", 2);
//		Calculate_Capacity_And_Ratio_V_over_C();
//
//		MainModual.WriteLog(0, "Step 2.10: Calculate Signal Delay", 2);
//		Calculate_Signal_Delay(nodeID);
//
//		MainModual.WriteLog(0, "Step 2.11: Judge LOS", 2);
//		Judge_Signal_LOS(nodeID);
//	}
//
//
//	void AddMovementStructure(int link_seq_no, string movement_str, int lanes, double signal_link_capacity, double saturation_flow_rate, int green_time,
//		int stage_no, int cycle_length, int priority, int sharedLanes, string linkID)
//	{
//		enum movement_Index mi = movement_str_to_index_map[movement_str];
//		enum direction di = movement_str_to_direction_map[movement_str];
//
//		movement_Array[mi].Enable = true;
//		movement_Array[mi].LinkSeqNo = link_seq_no;
//		//movement_Array[mi].Volume_per_lane = volume;//0
//		//movement_Array[mi].StageNo = stage(movement_to_Stage_Array[mi]);
//		movement_Array[mi].GroupNo = movement_Index_to_Group_Map[mi];
//		movement_Array[mi].DirectionNo = di;
//		movement_Array[mi].Left_Turn_Treatment = no_business;
//
//
//		movement_Array[mi].Lanes = lanes;
//		movement_Array[mi].signal_link_capacity = signal_link_capacity;
//		movement_Array[mi].saturation_flow_rate = saturation_flow_rate;
//		movement_Array[mi].green_time = green_time;
//		movement_Array[mi].pre_set_stage_no = stage_no;
//		movement_Array[mi].cycle_length = cycle_length;
//
//		//movement_Array[mi].priority = priority;
//		movement_Array[mi].SharedLanes = sharedLanes;
//		movement_Array[mi].linkID = linkID;
//		movement_Array[mi].Volume_per_lane = 500;
//	}
//
//	void AddMovementVolume_from_LinkPerformance()
//	{
//		CCSVParser_signal parser_link_Performance;
//		if (parser_link_Performance.OpenCSVFile("link_performance.csv", true))
//		{
//			while (parser_link_Performance.ReadRecord())  // if this line contains [] mark, then we will also read field headers.
//			{
//				string linkID;
//				parser_link_Performance.GetValueByFieldName("link_id", linkID);
//				for (size_t m = 1; m <= movement_Range; m++)
//				{
//					if (movement_Array[m].Enable == true && linkID == movement_Array[m].linkID)
//					{
//						int volume;
//						parser_link_Performance.GetValueByFieldName("D_per_hour_per_lane", volume);
//						movement_Array[m].Volume_per_lane = volume;
//					}
//				}
//			}
//		}
//		parser_link_Performance.CloseCSVFile();
//	}
//
//	void Set_Left_Turn_Treatment()
//	{
//		//20200808 Add a completeness check. Assign high privileges if only left.
//		for (size_t m = 1; m <= movement_Range; m++)
//		{
//			int final_decision;
//			if (movement_Array[m].GroupNo == 1)
//			{
//				final_decision = 0;
//				//(1)	Left-turn Lane Check
//				if (movement_Array[m].Lanes > 1)
//				{
//					final_decision = 1;
//				}
//				//(2)	Minimum Volume_per_lane Check
//				if (movement_Array[m].Volume_per_lane >= 240)
//				{
//					final_decision = 1;
//				}
//				//(3)	Opposing Through Lanes Check
//				int op_Movement_Index = left_Movement_Opposing_Index_Map[movement_Index(m)];
//				if (movement_Array[op_Movement_Index].Lanes >= 4)
//				{
//					final_decision = 1;
//				}
//				//(4)	Opposing Traffic Speed Check
//
//				//(5)	Minimum Cross-Product Check
//				int co_Movement_Index = left_Movement_Opposing_Index_Map[movement_Index(m)];
//				if (movement_Array[co_Movement_Index].Lanes > 1)
//				{
//					if (movement_Array[co_Movement_Index].Volume_per_lane * movement_Array[m].Volume_per_lane >= 100000)
//					{
//						final_decision = 1;
//					}
//				}
//				else
//				{
//
//					if (movement_Array[co_Movement_Index].Volume_per_lane * movement_Array[m].Volume_per_lane >= 50000)
//					{
//						final_decision = 1;
//					}
//				}
//				//(6) if there is no T movement, then the left movement should be protected.
//				if (movement_Array[left_Movement_Counterpart_Index_Map[movement_Index(m)]].Enable==false)
//				{
//					final_decision = 1;
//				}
//			}
//			else
//			{
//				final_decision = -1;
//			}
//			movement_Array[m].Left_Turn_Treatment = left_Turn_Treatment(final_decision);
//
//			if (final_decision != -1)
//			{
//				g_info_String = "Left Turn Treatment of Movement_";
//				g_info_String.append(movement_str_array[m]);
//				g_info_String.append(": ");
//				g_info_String.append(left_Turn_Treatment_index_to_str_map[left_Turn_Treatment(final_decision)]);
//				MainModual.WriteLog(1, g_info_String, 3);
//			}
//
//		}
//	}
//
//	void Set_StageNo_for_Movements()
//	{
//		stage_Range = 2;
//		int pre_set_num = 0;
//		for (int i = 0; i < movementSize; i++)
//		{
//			if (movement_Array[i].Enable = true && movement_Array[i].pre_set_stage_no > 0)
//			{
//				pre_set_num++;
//			}
//		}
//		if (pre_set_num > 0)
//		{
//			for (int i = 0; i < movementSize; i++)
//			{
//				if (movement_Array[i].Enable = true && movement_Array[i].pre_set_stage_no > 0)
//				{
//					movement_Array[i].StageNo_in_Order.push_back(stage(movement_Array[i].pre_set_stage_no));
//				}
//			}
//			return;
//
//		}
//
//
//		//Determining the main direction
//		int east_And_West_Volume = movement_Array[EBL].Volume_per_lane + movement_Array[EBT].Volume_per_lane + movement_Array[EBR].Volume_per_lane +
//			movement_Array[WBL].Volume_per_lane + movement_Array[WBT].Volume_per_lane + movement_Array[WBR].Volume_per_lane;
//
//		int north_And_South_Volume = movement_Array[NBL].Volume_per_lane + movement_Array[NBT].Volume_per_lane + movement_Array[NBR].Volume_per_lane +
//			movement_Array[SBL].Volume_per_lane + movement_Array[SBT].Volume_per_lane + movement_Array[SBR].Volume_per_lane;
//
//
//		bool east_And_West_Flag = false;
//		bool north_And_South_Flag = false;
//
//		if (movement_Array[EBL].Left_Turn_Treatment == prot)
//		{
//			stage_Range++;
//			east_And_West_Flag = true;
//		}
//		else if (movement_Array[WBL].Left_Turn_Treatment == prot)
//		{
//			stage_Range++;
//			east_And_West_Flag = true;
//
//		}
//
//		if (movement_Array[NBL].Left_Turn_Treatment == prot)
//		{
//			stage_Range++;
//			north_And_South_Flag = true;
//		}
//		else if (movement_Array[SBL].Left_Turn_Treatment == prot)
//		{
//			stage_Range++;
//			north_And_South_Flag = true;
//
//		}
//
//		enum movement_Index firstL;
//		enum movement_Index firstT;
//		enum movement_Index firstR;
//		enum movement_Index secondL;
//		enum movement_Index secondT;
//		enum movement_Index secondR;
//		enum movement_Index thridL;
//		enum movement_Index thridT;
//		enum movement_Index thridR;
//		enum movement_Index fouthL;
//		enum movement_Index fouthT;
//		enum movement_Index fouthR;
//		if (east_And_West_Volume >= north_And_South_Volume)
//		{
//			//east and west first
//			firstL = EBL;
//			firstT = EBT;
//			firstR = EBR;
//			secondL = WBL;
//			secondT = WBT;
//			secondR = WBR;
//			thridL = NBL;
//			thridT = NBT;
//			thridR = NBR;
//			fouthL = SBL;
//			fouthT = SBT;
//			fouthR = SBR;
//			MainModual.WriteLog(0, "Main Approaches: E & W", 3);
//		}
//		else
//		{
//
//			firstL = NBL;
//			firstT = NBT;
//			firstR = NBR;
//			secondL = SBL;
//			secondT = SBT;
//			secondR = SBR;
//			thridL = EBL;
//			thridT = EBT;
//			thridR = EBR;
//			fouthL = WBL;
//			fouthT = WBT;
//			fouthR = WBR;
//			MainModual.WriteLog(0, "Main Approaches: N & S", 3);
//
//		}
//
//
//		if (east_And_West_Flag)
//		{
//			movement_Array[firstL].StageNo_in_Order.push_back(stage1);
//			movement_Array[firstT].StageNo_in_Order.push_back(stage2);
//			movement_Array[firstR].StageNo_in_Order.push_back(stage2);
//			movement_Array[secondL].StageNo_in_Order.push_back(stage1);
//			movement_Array[secondT].StageNo_in_Order.push_back(stage2);
//			movement_Array[secondR].StageNo_in_Order.push_back(stage2);
//			if (north_And_South_Flag)
//			{
//				movement_Array[thridL].StageNo_in_Order.push_back(stage3);
//				movement_Array[thridT].StageNo_in_Order.push_back(stage4);
//				movement_Array[thridR].StageNo_in_Order.push_back(stage4);
//				movement_Array[fouthL].StageNo_in_Order.push_back(stage3);
//				movement_Array[fouthT].StageNo_in_Order.push_back(stage4);
//				movement_Array[fouthR].StageNo_in_Order.push_back(stage4);
//			}
//			else
//			{
//				movement_Array[thridL].StageNo_in_Order.push_back(stage3);
//				movement_Array[thridT].StageNo_in_Order.push_back(stage3);
//				movement_Array[thridR].StageNo_in_Order.push_back(stage3);
//				movement_Array[fouthL].StageNo_in_Order.push_back(stage3);
//				movement_Array[fouthT].StageNo_in_Order.push_back(stage3);
//				movement_Array[fouthR].StageNo_in_Order.push_back(stage3);
//			}
//
//		}
//		else
//		{
//			movement_Array[firstL].StageNo_in_Order.push_back(stage1);
//			movement_Array[firstT].StageNo_in_Order.push_back(stage1);
//			movement_Array[firstR].StageNo_in_Order.push_back(stage1);
//			movement_Array[secondL].StageNo_in_Order.push_back(stage1);
//			movement_Array[secondT].StageNo_in_Order.push_back(stage1);
//			movement_Array[secondR].StageNo_in_Order.push_back(stage1);
//			if (north_And_South_Flag)
//			{
//				movement_Array[thridL].StageNo_in_Order.push_back(stage2);
//				movement_Array[thridT].StageNo_in_Order.push_back(stage3);
//				movement_Array[thridR].StageNo_in_Order.push_back(stage3);
//				movement_Array[fouthL].StageNo_in_Order.push_back(stage2);
//				movement_Array[fouthT].StageNo_in_Order.push_back(stage3);
//				movement_Array[fouthR].StageNo_in_Order.push_back(stage3);
//			}
//			else
//			{
//				movement_Array[thridL].StageNo_in_Order.push_back(stage2);
//				movement_Array[thridT].StageNo_in_Order.push_back(stage2);
//				movement_Array[thridR].StageNo_in_Order.push_back(stage2);
//				movement_Array[fouthL].StageNo_in_Order.push_back(stage2);
//				movement_Array[fouthT].StageNo_in_Order.push_back(stage2);
//				movement_Array[fouthR].StageNo_in_Order.push_back(stage2);
//			}
//		}
//
//		//20200808 check the enable property, delete stage from enable=false movement
//		//20200812 modified this part by using vector StageNo_in_Order
//		int initialStage_Range = stage_Range;
//		int criticalStage = -1;
//
//		for (size_t s = initialStage_Range; s >= 1; s--)
//		{
//			int checkNumber = 0;
//			if (criticalStage != -1)
//			{
//				for (size_t m = 1; m <= movement_Range; m++)
//				{
//					for (size_t so = 0; so < movement_Array[m].StageNo_in_Order.size(); so++)
//					{
//						if (movement_Array[m].StageNo_in_Order[so] > criticalStage && movement_Array[m].Enable == true)
//						{
//							movement_Array[m].StageNo_in_Order[so] = stage(movement_Array[m].StageNo_in_Order[so] - 1);
//						}
//					}
//
//				}
//			}
//
//
//			for (size_t m = 1; m <= movement_Range; m++)
//			{
//				for (size_t so = 0; so < movement_Array[m].StageNo_in_Order.size(); so++)
//				{
//					if (movement_Array[m].Enable == true && movement_Array[m].StageNo_in_Order[so] == s)
//					{
//						checkNumber++;
//					}
//				}
//			}
//
//			if (checkNumber == 0)
//			{
//				stage_Range--;
//				criticalStage = s;
//			}
//			else
//			{
//				criticalStage = -1;
//			}
//		}
//		for (size_t m = 1; m <= movement_Range; m++)
//		{
//			if (movement_Array[m].Enable == false)
//			{
//				movement_Array[m].StageNo_in_Order[0] = stage(- 1);
//			}
//		}
//
//		//TODO: we can scan stages row-wise and mark the property of each stage, left-protected for example.
//
//		//for (size_t s = 0; s < stage_Range; s++)
//		//{
//
//		//}
//
//		//20200812 add right-turn treatment for movements 
//		for (size_t m = 1; m <= movement_Range; m++)
//		{
//			if (movement_Array[m].Enable==true && movement_Array[m].GroupNo==1)
//			{
//				if (movement_Array[left_Movement_Counterpart_Right_Trun_Index_Map[movement_Index(m)]].Enable == true)
//				{
//					movement_Array[left_Movement_Counterpart_Right_Trun_Index_Map[movement_Index(m)]].StageNo_in_Order.insert(movement_Array[left_Movement_Counterpart_Right_Trun_Index_Map[movement_Index(m)]].StageNo_in_Order.begin(), movement_Array[m].StageNo_in_Order[0]);
//				}
//			}
//
//		}
//
//
//
//		//movement_Array[EBL].StageNo = stage1;
//		//movement_Array[EBT].StageNo = stage2;
//		//movement_Array[EBR].StageNo = stage2;
//		//movement_Array[WBL].StageNo = stage1;
//		//movement_Array[WBT].StageNo = stage2;
//		//movement_Array[WBR].StageNo = stage2;
//		//movement_Array[NBL].StageNo = stage3;
//		//movement_Array[NBT].StageNo = stage3;
//		//movement_Array[NBR].StageNo = stage3;
//		//movement_Array[SBL].StageNo = stage3;
//		//movement_Array[SBT].StageNo = stage3;
//		//movement_Array[SBR].StageNo = stage3;
//
//		g_info_String = "Number of Stages: ";
//		g_info_String.append(to_string(stage_Range));
//		MainModual.WriteLog(1, g_info_String, 3);
//
//
//		for (size_t m = 1; m <= movement_Range; m++)
//		{
//			if (movement_Array[m].Enable == true)
//			{
//				for (size_t so = 0; so < movement_Array[m].StageNo_in_Order.size(); so++)
//				{
//					g_info_String = "StageNo of Movement_";
//					g_info_String.append(movement_str_array[m]);
//					g_info_String.append(": ");
//					g_info_String.append("Stage_");
//					g_info_String.append(to_string(movement_Array[m].StageNo_in_Order[so]));
//					MainModual.WriteLog(1, g_info_String, 4);
//				}
//			}
//
//		}
//	}
//
//	void Set_Saturation_Flow_Rate_Matrix()
//	{
//		//movement_Array left turn movement
//		// we need to use the saturation flow rate values based on protected and permitted
//
//		for (size_t m = 1; m <= movement_Range; m++)
//		{
//			if (movement_Array[m].Enable == false)
//			{
//				continue;
//			}
//			for (size_t so = 0; so < movement_Array[m].StageNo_in_Order.size(); so++)
//			{
//				if (movement_Array[m].Left_Turn_Treatment == prot)
//				{
//					saturation_Flow_Rate_Matrix[movement_Array[m].StageNo_in_Order[so]][m] = 1530 * movement_Array[m].Lanes * PHF;
//				}
//				else if (movement_Array[m].Left_Turn_Treatment == perm)
//				{
//					int op_Movement_Index = left_Movement_Opposing_Index_Map[movement_Index(m)];
//					int op_volume = movement_Array[op_Movement_Index].Volume_per_lane;
//					saturation_Flow_Rate_Matrix[movement_Array[m].StageNo_in_Order[so]][m] = f_1 * f_2 * op_volume * (exp(-op_volume * 4.5 / 3600)) / (1 - exp(-op_volume * 2.5 / 3600));
//				}
//				else
//				{
//					saturation_Flow_Rate_Matrix[movement_Array[m].StageNo_in_Order[so]][m] = 1530 * movement_Array[m].Lanes;//temp!!
//				}
//				g_info_String = "Saturation Flow Rate of Movement_";
//				g_info_String.append(movement_str_array[m]);
//				g_info_String.append(" and ");
//				g_info_String.append("Stage_");
//				g_info_String.append(to_string(movement_Array[m].StageNo_in_Order[so]));
//				g_info_String.append(": ");
//				g_info_String.append(to_string(saturation_Flow_Rate_Matrix[movement_Array[m].StageNo_in_Order[so]][m]));
//				MainModual.WriteLog(1, g_info_String, 3);
//			}
//		}
//
//		//saturation_Flow_Rate_Matrix[1][EBL] = 1750;
//		//saturation_Flow_Rate_Matrix[1][WBL] = 1750;
//		//saturation_Flow_Rate_Matrix[2][EBT] = 3400;
//		//saturation_Flow_Rate_Matrix[2][EBR] = 3400;
//		//saturation_Flow_Rate_Matrix[2][WBT] = 3400;
//		//saturation_Flow_Rate_Matrix[2][WBR] = 3400;
//		//saturation_Flow_Rate_Matrix[3][NBL] = 475;
//		//saturation_Flow_Rate_Matrix[3][NBT] = 1800;
//		//saturation_Flow_Rate_Matrix[3][NBR] = 1800;
//		//saturation_Flow_Rate_Matrix[3][SBL] = 450;
//		//saturation_Flow_Rate_Matrix[3][SBT] = 1800;
//		//saturation_Flow_Rate_Matrix[3][SBR] = 1800;
//	}
//
//	void Calculate_Flow_of_Ratio_Max()
//	{
//		//y_Stage_Movement_Matrix
//		//y_StageMax
//		for (size_t s = 1; s <= stage_Range; s++)
//		{
//			y_Max_Stage_Array[s] = 0;
//
//			for (size_t m = 1; m <= movement_Range; m++)
//			{
//				for (size_t so = 0; so < movement_Array[m].StageNo_in_Order.size(); so++)
//				{
//					if (saturation_Flow_Rate_Matrix[s][m] != 0 && movement_Array[m].Enable && movement_Array[m].StageNo_in_Order[so] == s)
//					{
//						y_Stage_Movement_Matrix[s][m] = double(movement_Array[m].Volume_per_lane) / double(saturation_Flow_Rate_Matrix[s][m]);
//
//						//double stage_Direction_Candidates_Matrix[stageSize][directionSize][groupSize]
//						stage_Direction_Candidates_Matrix[s][movement_Array[m].DirectionNo][movement_Array[m].GroupNo] += y_Stage_Movement_Matrix[s][m];
//
//						// we tally the movement matrix from this direction and this group number, so we can distingush movements belonging to different directions 
//						if (stage_Direction_Candidates_Matrix[s][movement_Array[m].DirectionNo][movement_Array[m].GroupNo] >= y_Max_Stage_Array[s])
//						{
//							y_Max_Stage_Array[s] = stage_Direction_Candidates_Matrix[s][movement_Array[m].DirectionNo][movement_Array[m].GroupNo];
//						}
//					}
//				}
//			}
//		}
//
//		y_StageMax = 0;
//		for (size_t i = 1; i <= stage_Range; i++)
//		{
//			y_StageMax += y_Max_Stage_Array[i];
//		}
//
//	}
//
//	void Calculate_Total_Cycle_Lost_Time()
//	{
//		l = t_L * stage_Range;
//	}
//
//	void Calculate_the_Minimum_And_Optimal_Cycle_Length()
//	{
//		c_Min = max(60, (l - x_c_Input) / (x_c_Input - y_StageMax));
//		c_Optimal= max(60, (1.5*l +5) / (1 - y_StageMax));
//	}
//
//	void Calculate_the_x_c_Output()
//	{
//		x_c_output = (y_StageMax * c_Min) / (c_Min - l);
//	}
//
//	void Calculate_Green_Time_for_Stages()
//	{
//		for (size_t s = 1; s <= stage_Range; s++)
//		{
////			green_Time_Stage_Array[i] = y_Max_Stage_Array[i] * c_Min / x_c_output;
//			green_Time_Stage_Array[s] = max(minGreenTime, y_Max_Stage_Array[s] * c_Min / y_StageMax);
//			effective_Green_Time_Stage_Array[s] = green_Time_Stage_Array[s] - t_L + t_Yellow + t_AR;
//			ratio_of_Effective_Green_Time_to_Cycle_Length_Array[s] = effective_Green_Time_Stage_Array[s] / c_Min;
//		}
//	}
//
//	void Printing_Green_Time_for_Stages()//output Effective green time 
//	{
//		cumulative_Green_Start_Time_Stage_Array[1] = 0;
//		cumulative_Green_End_Time_Stage_Array[1] = green_Time_Stage_Array[1];
//		MainModual.WriteLog(1, MainModual.Output_Combine_Message_And_Value_Double("Green Time of Stage 1", green_Time_Stage_Array[1]),3);
//		MainModual.WriteLog(1, MainModual.Output_Combine_Message_And_Value_Double("Start Green Time of Stage 1", cumulative_Green_Start_Time_Stage_Array[1]), 4);
//		MainModual.WriteLog(1, MainModual.Output_Combine_Message_And_Value_Double("End Green Time of Stage 1", cumulative_Green_End_Time_Stage_Array[1]), 4);
//
//
//		for (size_t i = 2; i <= stage_Range; i++)
//		{
//			g_info_String = "Green Time of Stage ";
//			g_info_String.append(to_string(i));
//			MainModual.WriteLog(1, MainModual.Output_Combine_Message_And_Value_Double(g_info_String, green_Time_Stage_Array[i]),3);
//			cumulative_Green_Start_Time_Stage_Array[i] = cumulative_Green_End_Time_Stage_Array[i - 1];
//			cumulative_Green_End_Time_Stage_Array[i] = cumulative_Green_Start_Time_Stage_Array[i] + green_Time_Stage_Array[i];
//
//			g_info_String = "Start Green Time of Stage ";
//			g_info_String.append(to_string(i));
//			MainModual.WriteLog(1, MainModual.Output_Combine_Message_And_Value_Double(g_info_String, cumulative_Green_Start_Time_Stage_Array[i]), 4);
//
//			g_info_String = "End Green Time of Stage ";
//			g_info_String.append(to_string(i));
//			MainModual.WriteLog(1, MainModual.Output_Combine_Message_And_Value_Double(g_info_String, cumulative_Green_End_Time_Stage_Array[i]), 4);
//		}
//
//
//
//		cumulative_Effective_Green_Start_Time_Stage_Array[1] = 0;
//		cumulative_Effective_Green_End_Time_Stage_Array[1] = effective_Green_Time_Stage_Array[1];
//
//		MainModual.WriteLog(1, MainModual.Output_Combine_Message_And_Value_Double("Effective Green Time of Stage 1", effective_Green_Time_Stage_Array[1]), 3);
//		MainModual.WriteLog(1, MainModual.Output_Combine_Message_And_Value_Double("Start Effective Green Time of Stage 1", cumulative_Effective_Green_Start_Time_Stage_Array[1]), 4);
//		MainModual.WriteLog(1, MainModual.Output_Combine_Message_And_Value_Double("End Effective Green Time of Stage 1", cumulative_Effective_Green_End_Time_Stage_Array[1]), 4);
//
//
//		for (size_t i = 2; i <= stage_Range; i++)
//		{
//			g_info_String = "Effective Green Time of Stage ";
//			g_info_String.append(to_string(i));
//			MainModual.WriteLog(1, MainModual.Output_Combine_Message_And_Value_Double(g_info_String, effective_Green_Time_Stage_Array[i]), 3);
//			cumulative_Effective_Green_Start_Time_Stage_Array[i] = cumulative_Effective_Green_End_Time_Stage_Array[i - 1];
//			cumulative_Effective_Green_End_Time_Stage_Array[i] = cumulative_Effective_Green_Start_Time_Stage_Array[i] + effective_Green_Time_Stage_Array[i];
//
//			g_info_String = "Start Effective Green Time of Stage ";
//			g_info_String.append(to_string(i));
//			MainModual.WriteLog(1, MainModual.Output_Combine_Message_And_Value_Double(g_info_String, cumulative_Effective_Green_Start_Time_Stage_Array[i]), 4);
//
//			g_info_String = "End Effective Green Time of Stage ";
//			g_info_String.append(to_string(i));
//			MainModual.WriteLog(1, MainModual.Output_Combine_Message_And_Value_Double(g_info_String, cumulative_Effective_Green_End_Time_Stage_Array[i]), 4);
//		}
//	}
//
//	void Calculate_Capacity_And_Ratio_V_over_C()
//	{
//		for (size_t s = 1; s <= stage_Range; s++)
//		{
//			for (size_t m = 1; m <= movement_Range; m++)
//			{
//				for (size_t so = 0; so < movement_Array[m].StageNo_in_Order.size(); so++)
//				{
//					if (saturation_Flow_Rate_Matrix[s][m] != 0 && movement_Array[m].Enable && movement_Array[m].StageNo_in_Order[so] == s)
//					{
//						capacity_by_Stage_and_Movement_Matrix[s][m] = saturation_Flow_Rate_Matrix[s][m] * ratio_of_Effective_Green_Time_to_Cycle_Length_Array[s];
//						g_info_String = "a. Capacity of Stage_";
//						g_info_String.append(to_string(s));
//						g_info_String.append(" and Movement_");
//						g_info_String.append(movement_str_array[m]);
//						g_info_String.append(": ");
//						g_info_String.append(to_string(capacity_by_Stage_and_Movement_Matrix[s][m]));
//						MainModual.WriteLog(1, g_info_String, 3);
//						v_over_C_by_Stage_and_Movement_Matrix[s][m] = movement_Array[m].Volume_per_lane / capacity_by_Stage_and_Movement_Matrix[s][m];
//						g_info_String = "b. V/C of Stage_";
//						g_info_String.append(to_string(s));
//						g_info_String.append(" and Movement_");
//						g_info_String.append(movement_str_array[m]);
//						g_info_String.append(": ");
//						g_info_String.append(to_string(v_over_C_by_Stage_and_Movement_Matrix[s][m]));
//						MainModual.WriteLog(1, g_info_String, 3);
//					}
//				}
//			}
//		}
//	}
//
//	void Calculate_Signal_Delay(int nodeID)
//	{
//		for (size_t s = 1; s <= stage_Range; s++)
//		{
//			for (size_t m = 1; m <= movement_Range; m++)
//			{
//				for (size_t so = 0; so < movement_Array[m].StageNo_in_Order.size(); so++)
//				{
//					if (saturation_Flow_Rate_Matrix[s][m] != 0 && movement_Array[m].Enable && movement_Array[m].StageNo_in_Order[so] == s)
//					{
//						average_Uniform_Delay_Matrix[s][m] = 
//							(0.5 * capacity_by_Stage_and_Movement_Matrix[s][m] * pow((1 - ratio_of_Effective_Green_Time_to_Cycle_Length_Array[s]), 2)) 
//							/ (1 - v_over_C_by_Stage_and_Movement_Matrix[s][m] * ratio_of_Effective_Green_Time_to_Cycle_Length_Array[s]);
//
//						//(2)	Average Incremental Delay          next time 
//						g_info_String = "Average Uniform Delay of Stage_";
//						g_info_String.append(to_string(s));
//						g_info_String.append(" and Movement_");
//						g_info_String.append(movement_str_array[m]);
//						g_info_String.append(": ");
//						g_info_String.append(to_string(average_Uniform_Delay_Matrix[s][m]));
//						MainModual.WriteLog(1, g_info_String, 3);
//
//
//						approach_Total_Delay_Array[movement_Array[m].DirectionNo] += movement_Array[m].Volume_per_lane * average_Uniform_Delay_Matrix[s][m];
//						approach_Total_Volume_Array[movement_Array[m].DirectionNo] += movement_Array[m].Volume_per_lane;
//					}
//				}
//			}
//		}
//		for (size_t d = 1; d <= direction_Range; d++)
//		{
//			if (approach_Total_Volume_Array[d] == 0)
//			{
//				continue;
//			}
//			approach_Average_Delay_Array[d] = approach_Total_Delay_Array[d] / approach_Total_Volume_Array[d];
//			g_info_String = "Average Delay of Approach_";
//			g_info_String.append(direction_index_to_str_map[direction(d)]);
//			g_info_String.append(": ");
//			g_info_String.append(to_string(approach_Average_Delay_Array[d]));
//			MainModual.WriteLog(1, g_info_String, 3);
//
//
//
//			intersection_Total_Delay += approach_Average_Delay_Array[d] * approach_Total_Volume_Array[d];
//			intersection_Total_Volume += approach_Total_Volume_Array[d];
//		}
//		intersection_Average_Delay = intersection_Total_Delay / intersection_Total_Volume;
//
//		g_info_String = "Total Delay of Intersection_NodeID_";
//		g_info_String.append(to_string(nodeID));
//		g_info_String.append(": ");
//		g_info_String.append(to_string(intersection_Total_Delay));
//		MainModual.WriteLog(1, g_info_String, 3);
//
//
//
//		g_info_String = "Average Delay of Intersection_NodeID_";
//		g_info_String.append(to_string(nodeID));
//		g_info_String.append(": ");
//		g_info_String.append(to_string(intersection_Average_Delay));
//		MainModual.WriteLog(1, g_info_String, 3);
//	}
//
//	void Judge_Signal_LOS(int nodeID)
//	{
//		if (intersection_Total_Delay <= 10)
//		{
//			LOS = "A";
//		}
//		else if (intersection_Total_Delay <= 20)
//		{
//			LOS = "B";
//		}
//		else if (intersection_Total_Delay <= 35)
//		{
//			LOS = "C";
//		}
//		else if (intersection_Total_Delay <= 55)
//		{
//			LOS = "D";
//		}
//		else if (intersection_Total_Delay <= 80)
//		{
//			LOS = "E";
//		}
//		else
//		{
//			LOS = "F";
//		}
//		g_info_String = "LOS of Intersection_NodeID_";
//		g_info_String.append(to_string(nodeID));
//		g_info_String.append(": ");
//		g_info_String.append(LOS);
//		MainModual.WriteLog(1, g_info_String, 3);
//	}
//
//	int signal_node_seq_no;  // sequence number 
//	int main_node_id;      //external node number 
//
//	std::vector<int> m_movement_link_seq_no_vector;
//
//};
//
//std::map<int, C_SignalNode> g_signal_node_map;  // first key is signal node id
//std::vector<C_SignalNode> g_signal_node_vector;
//std::vector<C_SignalLink> g_signal_link_vector;
//
//
//vector<string> split_sa(const string& s, const string& seperator) {
//	vector<string> result;
//	typedef string::size_type string_size;
//	string_size i = 0;
//
//	while (i != s.size()) {
//		int flag = 0;
//		while (i != s.size() && flag == 0) {
//			flag = 1;
//			for (string_size x = 0; x < seperator.size(); ++x)
//				if (s[i] == seperator[x]) {
//					++i;
//					flag = 0;
//					break;
//				}
//		}
//
//		flag = 0;
//		string_size j = i;
//		while (j != s.size() && flag == 0) {
//			for (string_size x = 0; x < seperator.size(); ++x)
//				if (s[j] == seperator[x]) {
//					flag = 1;
//					break;
//				}
//			if (flag == 0)
//				++j;
//		}
//		if (i != j) {
//			result.push_back(s.substr(i, j - i));
//			i = j;
//		}
//	}
//	return result;
//}
//
//
//vector<float> g_sa_time_parser(vector<string>& inputstring)
//{
//	vector<float> output_global_minute;
//
//	for (int k = 0; k < inputstring.size(); k++)
//	{
//		vector<string> sub_string = split_sa(inputstring[k], "_");
//
//		for (int i = 0; i < sub_string.size(); i++)
//		{
//			//HHMM
//			//012345
//			char hh1 = sub_string[i].at(0);
//			char hh2 = sub_string[i].at(1);
//			char mm1 = sub_string[i].at(2);
//			char mm2 = sub_string[i].at(3);
//
//			float hhf1 = ((float)hh1 - 48);
//			float hhf2 = ((float)hh2 - 48);
//			float mmf1 = ((float)mm1 - 48);
//			float mmf2 = ((float)mm2 - 48);
//
//			float hh = hhf1 * 10 * 60 + hhf2 * 60;
//			float mm = mmf1 * 10 + mmf2;
//			float global_mm_temp = hh + mm;
//			output_global_minute.push_back(global_mm_temp);
//		}
//	}
//
//	return output_global_minute;
//} // transform hhmm to minutes 
//
//
//void sa_ReadMesoData(C_SignalMainModual& MainModual)
//{
//	MainModual.g_number_of_nodes = 0;
//	MainModual.g_number_of_links = 0;  // initialize  the counter to 0
//
//	int internal_node_seq_no = 0;
//	// step X: read node file 
//
//	CCSVParser_signal parser_node;
//	if (parser_node.OpenCSVFile("node.csv", true))
//	{
//		while (parser_node.ReadRecord())  // if this line contains [] mark, then we will also read field headers.
//		{
//			int node_id;
//			if (parser_node.GetValueByFieldName("node_id", node_id) == false)
//				continue;
//
//			if (MainModual.g_internal_node_to_seq_no_map.find(node_id) != MainModual.g_internal_node_to_seq_no_map.end())
//			{
//				continue; //has been defined
//			}
//			MainModual.g_internal_node_to_seq_no_map[node_id] = internal_node_seq_no;
//			C_SignalNode node;  // create a node object 
//			node.node_id = node_id;
//			node.node_seq_no = internal_node_seq_no;
//			internal_node_seq_no++;
//			g_signal_node_vector.push_back(node);  // push it to the global node vector
//			MainModual.g_number_of_nodes++;
//		}
//		g_info_String = "Number of Nodes = ";
//		g_info_String.append(to_string(MainModual.g_number_of_nodes));
//		MainModual.WriteLog(0, g_info_String, 2);
//		parser_node.CloseCSVFile();
//	}
//
//	// step 0: read link file 
//
//	CCSVParser_signal parser_link;
//	if (parser_link.OpenCSVFile("link.csv", true))
//	{
//		while (parser_link.ReadRecord())  // if this line contains [] mark, then we will also read field headers.
//		{
//			int from_node_id;
//			int to_node_id;
//			if (parser_link.GetValueByFieldName("from_node_id", from_node_id) == false)
//				continue;
//			if (parser_link.GetValueByFieldName("to_node_id", to_node_id) == false)
//				continue;
//			string linkID;
//			if (parser_link.GetValueByFieldName("link_id", linkID) == false)
//				continue;
//
//			// add the to node id into the outbound (adjacent) node list
//
//			if (MainModual.g_internal_node_to_seq_no_map.find(from_node_id) == MainModual.g_internal_node_to_seq_no_map.end())
//			{
//				cout << "Error: from_node_id " << from_node_id << " in file road_link.csv is not defined in node.csv." << endl;
//				continue; //has not been defined
//			}
//			if (MainModual.g_internal_node_to_seq_no_map.find(to_node_id) == MainModual.g_internal_node_to_seq_no_map.end())
//			{
//				cout << "Error: to_node_id " << to_node_id << " in file road_link.csv is not defined in node.csv." << endl;
//				continue; //has not been defined
//			}
//
//			if (MainModual.g_road_link_id_map.find(linkID) != MainModual.g_road_link_id_map.end())
//			{
//				cout << "Error: road_link_id " << linkID.c_str() << " has been defined more than once. Please check road_link.csv." << endl;
//				continue; //more than once.
//			}
//
//			int internal_from_node_seq_no = MainModual.g_internal_node_to_seq_no_map[from_node_id];  // map external node number to internal node seq no. 
//			int internal_to_node_seq_no = MainModual.g_internal_node_to_seq_no_map[to_node_id];
//
//
//			C_SignalLink link;  // create a link object 
//
//			link.from_node_seq_no = internal_from_node_seq_no;
//			link.to_node_seq_no = internal_to_node_seq_no;
//			link.link_seq_no = MainModual.g_number_of_links;
//			link.to_node_seq_no = internal_to_node_seq_no;
//			link.link_id = linkID;
//
//			MainModual.g_road_link_id_map[link.link_id] = 1;
//
//			parser_link.GetValueByFieldName("facility_type", link.link_type, true, false);
//			parser_link.GetValueByFieldName("link_type", link.link_type);
//
//
//			float length = 1.0; // km or mile
//			float free_speed = 1.0;
//
//			float signal_link_capacity = 1800;
//			parser_link.GetValueByFieldName("length", length);
//			parser_link.GetValueByFieldName("free_speed", free_speed);
//			free_speed = max(0.1, free_speed);
//
//			int number_of_lanes = 1;
//			parser_link.GetValueByFieldName("lanes", number_of_lanes);
//			parser_link.GetValueByFieldName("signal_link_capacity", signal_link_capacity);
//
//			float default_cap = 1000;
//			float default_BaseTT = 1;
//
//			link.number_of_lanes = number_of_lanes;
//			link.lane_capacity = signal_link_capacity;
//
//
//			link.length = length;
//			link.free_flow_travel_time_in_min = length / free_speed * 60;
//
//
//			g_signal_node_vector[internal_from_node_seq_no].m_outgoing_link_seq_no_vector.push_back(link.link_seq_no);  // add this link to the corresponding node as part of outgoing node/link
//			g_signal_node_vector[internal_to_node_seq_no].m_incoming_link_seq_no_vector.push_back(link.link_seq_no);  // add this link to the corresponding node as part of outgoing node/link
//
//			g_signal_node_vector[internal_from_node_seq_no].m_to_node_seq_no_vector.push_back(link.to_node_seq_no);  // add this link to the corresponding node as part of outgoing node/link
//			g_signal_node_vector[internal_from_node_seq_no].m_to_node_2_link_seq_no_map[link.to_node_seq_no] = link.link_seq_no;  // add this link to the corresponding node as part of outgoing node/link
//			g_signal_link_vector.push_back(link);
//
//			MainModual.g_number_of_links++;
//
//			// map link data to signal node map.
//
//			string movement_str;
//			parser_link.GetValueByFieldName("movement_str", movement_str);
//
//
//			if (movement_str.size() > 0)  // and valid
//			{
//				int main_node_id = -1;
//				parser_link.GetValueByFieldName("main_node_id", main_node_id);
//
//				int lanes = 0;
//				parser_link.GetValueByFieldName("lanes", lanes);
//
//
//				double saturation_flow_rate = 0;
//				//parser_link.GetValueByFieldName("saturation_flow_rate", saturation_flow_rate);
//
//				int green_time = 0;
//				parser_link.GetValueByFieldName("green_time", green_time);
//
//				int stage_no = -1;
//				parser_link.GetValueByFieldName("stage_no", stage_no);
//
//				int cycle_length = 0;
//				parser_link.GetValueByFieldName("cycle_length", cycle_length);
//
//				int priority = 0;
//				//parser_link.GetValueByFieldName("priority", priority);
//
//				int sharedLanes = 0;
//				//parser_link.GetValueByFieldName("sharedLanes", sharedLanes);
//
//				if (main_node_id >= 1)
//				{
//					g_info_String = "Main Node ID: ";
//					g_info_String.append(to_string(main_node_id));
//					g_info_String.append("| movement_str: ");
//					g_info_String.append(movement_str);
//					MainModual.WriteLog(0, g_info_String, 2);
//
//
//
//					g_signal_node_map[main_node_id].AddMovementStructure(link.link_seq_no, movement_str,
//						lanes, signal_link_capacity, saturation_flow_rate, green_time, stage_no, cycle_length, priority, sharedLanes, linkID);// we add parameters here
//				}
//			}
//		}
//	}
//		parser_link.CloseCSVFile();
//
//	// we now know the number of links
//		g_info_String = "Number of Links = ";
//		g_info_String.append(to_string(MainModual.g_number_of_links));
//		MainModual.WriteLog(0, g_info_String, 2);
//};
//
//
//std::vector<C_SignalNode> SignalAPI(int iteration_number, int MainModual_mode, int column_updating_iterations)
//{
//
//	// step 1: read input data of network / demand tables / Toll
//	sa_ReadMesoData(MainModual);
//	for (std::map<int, C_SignalNode>::iterator it = g_signal_node_map.begin(); it != g_signal_node_map.end(); ++it)
//	{
//		it->second.PerformQEM(it->first);
//
//		//output mainnode ID
//	}
//
//
//	for (std::map<int, C_SignalNode>::iterator it = g_signal_node_map.begin(); it != g_signal_node_map.end(); ++it)
//	{
//		C_SignalNode sn = it->second;
//
//		int cycle_time_in_sec = max(10, sn.c_Min);
//		int number_of_cycles = (MainModual.g_LoadingEndTimeInMin - MainModual.g_LoadingStartTimeInMin) * 60 / cycle_time_in_sec;  // unit: seconds;
//		int offset_in_sec = 0;
//		int g_loading_start_time_in_sec = MainModual.g_LoadingStartTimeInMin * 60 + offset_in_sec;
//		int ci = 0;
//
//		for (int m = 1; m < movementSize; m++)
//		{
//			if (sn.movement_Array[m].Enable)
//			{
//				for (size_t so = 0; so < sn.movement_Array[m].StageNo_in_Order.size(); so++)
//				{
//					int StageNo = sn.movement_Array[m].StageNo_in_Order[so];
//					// we should also consider offset.
//					int global_start_time_in_sec = sn.cumulative_Green_Start_Time_Stage_Array[StageNo] + cycle_time_in_sec * ci + g_loading_start_time_in_sec;
//					int global_end_time_in_sec = sn.cumulative_Green_End_Time_Stage_Array[StageNo] + cycle_time_in_sec * ci + g_loading_start_time_in_sec;
//
//					//0300:30
//					int start_hour = global_start_time_in_sec / 3600;
//					int start_min = global_start_time_in_sec / 60 - start_hour * 60;
//					int start_sec = global_start_time_in_sec % 60;
//
//					int end_hour = global_end_time_in_sec / 3600;
//					int end_min = global_end_time_in_sec / 60 - end_hour * 60;
//					int end_sec = global_end_time_in_sec % 60;
//
//					int from_node_id = g_signal_node_vector[g_signal_link_vector[sn.movement_Array[m].LinkSeqNo].from_node_seq_no].node_id;
//					int to_node_id = g_signal_node_vector[g_signal_link_vector[sn.movement_Array[m].LinkSeqNo].to_node_seq_no].node_id;
//					//						float capacity = sn.green_Time_Stage_Array[StageNo] * sn.saturation_Flow_Rate_Matrix[StageNo][m] / 3600.0;
//					float capacity = sn.green_Time_Stage_Array[StageNo] * 1800.0 / 3600.0;
//					//float capacity = sn.capacity_by_Stage_and_Movement_Matrix[StageNo][m]/60;
//
//					float greenTime = sn.green_Time_Stage_Array[StageNo];
//					g_signal_link_vector[sn.movement_Array[m].LinkSeqNo].lane_capacity = capacity;
//					g_signal_link_vector[sn.movement_Array[m].LinkSeqNo].greenTime = greenTime;
//				}
//			}
//		}
//	}
//	MainModual.WriteLog(0, "-----------------------Finished----------------------- ", 1);
//	return g_signal_node_vector;// we return sa_node_vector
//}