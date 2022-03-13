/* Portions Copyright 2021 Xuesong Zhou, Peiheng Li, Cafer Avci
 *
 * If you help write or modify the code, please also list your names here.
 * The reason of having Copyright info here is to ensure all the modified version, as a whole, under the GPL
 * and further prevent a violation of the GPL.
 *
 * More about "How to use GNU licenses for your own software"
 * http://www.gnu.org/licenses/gpl-howto.html
 */

#include "teestream.h"
#include "config.h"
#include "utils.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <string>
#include <cstring>
#include <vector>
#include <cmath>

using std::endl;
using std::string;
using std::vector;
using std::ofstream;
using std::istringstream;
using std::ostringstream;

using std::asin;
using std::sin;
using std::cos;
using std::pow;
using std::sqrt;
using std::min;
using std::fmin;

void g_program_stop()
{
    dtalog.output() << "DTALite Program stops. Press any key to terminate. Thanks!" << endl;
    getchar();
    exit(0);
}

void g_program_exit()
{
    dtalog.output() << "DTALite Program completes. Thanks!" << endl;

    exit(0);
}


void fopen_ss(FILE** file, const char* fileName, const char* mode)
{
    *file = fopen(fileName, mode);
}

float g_read_float(FILE* f)
{
    if (feof(f) == 1)
        return -1;
    /*
        read a floating point number from the current pointer of the file,
        skip all spaces
     */
    char ch, buf[32];
    int i = 0;
    int flag = 1;

    /* returns -1 if end of file is reached */
    while (true)
    {
        ch = getc(f);
        if (ch == EOF || ch == '*' || ch == '$' || ch < -1 || ch >=255 )
            return -1;

        if (isdigit(ch))
            break;

        if (ch == '-')
            flag = -1;
        else
            flag = 1;
    }

    if (ch == EOF) return -1;
    while (isdigit(ch) || ch == '.') {
        buf[i++] = ch;
        ch = fgetc(f);

    }
    buf[i] = 0;

    /* atof function converts a character string (char *) into a doubleing
    pointer equivalent, and if the string is not a floting point number,
    a zero will be return.
    */

    return (float)(atof(buf) * flag);
}

//split the string by "_"
vector<string> split(const string &s, const string &seperator)
{
    vector<string> result;
    typedef string::size_type string_size;
    string_size i = 0;

    while (i != s.size()) {
        int flag = 0;
        while (i != s.size() && flag == 0) {
            flag = 1;
            for (string_size x = 0; x < seperator.size(); ++x)
                if (s[i] == seperator[x]) {
                    ++i;
                    flag = 0;
                    break;
                }
        }

        flag = 0;
        string_size j = i;
        while (j != s.size() && flag == 0) {
            for (string_size x = 0; x < seperator.size(); ++x)
                if (s[j] == seperator[x]) {
                    flag = 1;
                    break;
                }
            if (flag == 0)
                ++j;
        }
        if (i != j) {
            result.push_back(s.substr(i, j - i));
            i = j;
        }
    }

    return result;
}

vector<float> g_time_parser(string str)
{
    vector<float> output_global_minute;

    int string_lenghth = str.length();

    //ASSERT(string_lenghth < 100);

    const char* string_line = str.data(); //string to char*

    int char_length = strlen(string_line);

    char ch, buf_ddhhmm[32] = { 0 }, buf_SS[32] = { 0 }, buf_sss[32] = { 0 };
    char dd1, dd2, hh1, hh2, mm1, mm2, SS1, SS2, sss1, sss2, sss3;
    float ddf1, ddf2, hhf1, hhf2, mmf1, mmf2, SSf1, SSf2, sssf1, sssf2, sssf3;
    float global_minute = 0;
    float dd = 0, hh = 0, mm = 0, SS = 0, sss = 0;
    int i = 0;
    int buffer_i = 0, buffer_k = 0, buffer_j = 0;
    int num_of_colons = 0;

    //DDHHMM:SS:sss or HHMM:SS:sss

    while (i < char_length)
    {
        ch = string_line[i++];

        if (num_of_colons == 0 && ch != '_' && ch != ':') //input to buf_ddhhmm until we meet the colon
        {
            buf_ddhhmm[buffer_i++] = ch;
        }
        else if (num_of_colons == 1 && ch != ':') //start the Second "SS"
        {
            buf_SS[buffer_k++] = ch;
        }
        else if (num_of_colons == 2 && ch != ':') //start the Millisecond "sss"
        {
            buf_sss[buffer_j++] = ch;
        }

        if (ch == '_' || i == char_length) //start a new time string
        {
            if (buffer_i == 4) //"HHMM"
            {
                //HHMM, 0123
                hh1 = buf_ddhhmm[0]; //read each first
                hh2 = buf_ddhhmm[1];
                mm1 = buf_ddhhmm[2];
                mm2 = buf_ddhhmm[3];

                hhf1 = ((float)hh1 - 48); //convert a char to a float
                hhf2 = ((float)hh2 - 48);
                mmf1 = ((float)mm1 - 48);
                mmf2 = ((float)mm2 - 48);

                dd = 0;
                hh = hhf1 * 10 * 60 + hhf2 * 60;
                mm = mmf1 * 10 + mmf2;
            }
            else if (buffer_i == 6) //"DDHHMM"
            {
                //DDHHMM, 012345
                dd1 = buf_ddhhmm[0]; //read each first
                dd2 = buf_ddhhmm[1];
                hh1 = buf_ddhhmm[2];
                hh2 = buf_ddhhmm[3];
                mm1 = buf_ddhhmm[4];
                mm2 = buf_ddhhmm[5];

                ddf1 = ((float)dd1 - 48); //convert a char to a float
                ddf2 = ((float)dd2 - 48);
                hhf1 = ((float)hh1 - 48);
                hhf2 = ((float)hh2 - 48);
                mmf1 = ((float)mm1 - 48);
                mmf2 = ((float)mm2 - 48);

                dd = ddf1 * 10 * 24 * 60 + ddf2 * 24 * 60;
                hh = hhf1 * 10 * 60 + hhf2 * 60;
                mm = mmf1 * 10 + mmf2;
            }

            if (num_of_colons == 1 || num_of_colons == 2)
            {
                //SS, 01
                SS1 = buf_SS[0]; //read each first
                SS2 = buf_SS[1];

                SSf1 = ((float)SS1 - 48); //convert a char to a float
                SSf2 = ((float)SS2 - 48);

                SS = (SSf1 * 10 + SSf2) / 60;
            }

            if (num_of_colons == 2)
            {
                //sss, 012
                sss1 = buf_sss[0]; //read each first
                sss2 = buf_sss[1];
                sss3 = buf_sss[2];

                sssf1 = ((float)sss1 - 48); //convert a char to a float
                sssf2 = ((float)sss2 - 48);
                sssf3 = ((float)sss3 - 48);

                sss = (sssf1 * 100 + sssf2 * 10 + sssf3) / 1000;
            }

            global_minute = dd + hh + mm + SS + sss;

            output_global_minute.push_back(global_minute);

            //initialize the parameters
            buffer_i = 0;
            buffer_k = 0;
            buffer_j = 0;
            num_of_colons = 0;
        }

        if (ch == ':')
            num_of_colons += 1;
    }

    return output_global_minute;
}

float g_timestamp_parser(string str)
{
   float output_global_minute;

    int string_lenghth = str.length();

    //ASSERT(string_lenghth < 100);

    const char* string_line = str.data(); //string to char*

    int char_length = strlen(string_line);

    char ch, buf_ddhhmm[32] = { 0 }, buf_SS[32] = { 0 }, buf_sss[32] = { 0 };
    char dd1, dd2, hh1, hh2, mm1, mm2, SS1, SS2, sss1, sss2, sss3;
    float ddf1, ddf2, hhf1, hhf2, mmf1, mmf2, SSf1, SSf2, sssf1, sssf2, sssf3;
    float global_minute = 0;
    float dd = 0, hh = 0, mm = 0, SS = 0, sss = 0;
    int i = 1;  // skip T as the first letter
    int buffer_i = 0, buffer_k = 0, buffer_j = 0;
    int num_of_colons = 0;

    //DDHHMM:SS:sss or HHMM:SS:sss

    while (i < char_length)
    {
        ch = string_line[i++];

        if (num_of_colons == 0 && ch != '_' && ch != ':') //input to buf_ddhhmm until we meet the colon
        {
            buf_ddhhmm[buffer_i++] = ch;
        }
        else if (num_of_colons == 1 && ch != ':') //start the Second "SS"
        {
            buf_SS[buffer_k++] = ch;
        }
        else if (num_of_colons == 2 && ch != ':') //start the Millisecond "sss"
        {
            buf_sss[buffer_j++] = ch;
        }

        if (i == char_length) //start a new time string
        {
            if (buffer_i == 4) //"HHMM"
            {
                //HHMM, 0123
                hh1 = buf_ddhhmm[0]; //read each first
                hh2 = buf_ddhhmm[1];
                mm1 = buf_ddhhmm[2];
                mm2 = buf_ddhhmm[3];

                hhf1 = ((float)hh1 - 48); //convert a char to a float
                hhf2 = ((float)hh2 - 48);
                mmf1 = ((float)mm1 - 48);
                mmf2 = ((float)mm2 - 48);

                dd = 0;
                hh = hhf1 * 10 * 60 + hhf2 * 60;
                mm = mmf1 * 10 + mmf2;
            }
            else if (buffer_i == 6) //"DDHHMM"
            {
                //DDHHMM, 012345
                dd1 = buf_ddhhmm[0]; //read each first
                dd2 = buf_ddhhmm[1];
                hh1 = buf_ddhhmm[2];
                hh2 = buf_ddhhmm[3];
                mm1 = buf_ddhhmm[4];
                mm2 = buf_ddhhmm[5];

                ddf1 = ((float)dd1 - 48); //convert a char to a float
                ddf2 = ((float)dd2 - 48);
                hhf1 = ((float)hh1 - 48);
                hhf2 = ((float)hh2 - 48);
                mmf1 = ((float)mm1 - 48);
                mmf2 = ((float)mm2 - 48);

                dd = ddf1 * 10 * 24 * 60 + ddf2 * 24 * 60;
                hh = hhf1 * 10 * 60 + hhf2 * 60;
                mm = mmf1 * 10 + mmf2;
            }

            if (num_of_colons == 1 || num_of_colons == 2)
            {
                //SS, 01
                SS1 = buf_SS[0]; //read each first
                SS2 = buf_SS[1];

                SSf1 = ((float)SS1 - 48); //convert a char to a float
                SSf2 = ((float)SS2 - 48);

                SS = (SSf1 * 10 + SSf2) / 60;
            }

            if (num_of_colons == 2)
            {
                //sss, 012
                sss1 = buf_sss[0]; //read each first
                sss2 = buf_sss[1];
                sss3 = buf_sss[2];

                sssf1 = ((float)sss1 - 48); //convert a char to a float
                sssf2 = ((float)sss2 - 48);
                sssf3 = ((float)sss3 - 48);

                sss = (sssf1 * 100 + sssf2 * 10 + sssf3) / 1000;
            }

            global_minute = dd + hh + mm + SS + sss;

            output_global_minute = global_minute;

            //initialize the parameters
            buffer_i = 0;
            buffer_k = 0;
            buffer_j = 0;
            num_of_colons = 0;
        }

        if (ch == ':')
            num_of_colons += 1;
    }

    return output_global_minute;
}

string g_time_coding(float time_stamp)
{
    int hour = static_cast<int>(time_stamp / 60);
    int minute = static_cast<int>(time_stamp - hour * 60);
    int second = static_cast<int>((time_stamp - hour * 60 - minute) * 60 + 0.02);

    int sss = ((time_stamp - hour * 60 - minute) * 60 - second)*1000;

    //mm:ss.sss
    ostringstream strm;
    strm.fill('0');
    strm << std::setw(2) << hour << std::setw(2) << minute << ":" << std::setw(2) << second << "." << std::setw(3) << sss;

    return strm.str();
}


int g_ParserStringSequence(std::string str_input, vector<string>& vect)
{

    std::istringstream ss(str_input);
    std::string token;

    while (std::getline(ss, token, ';')) {
        vect.push_back(token);
    }

    return vect.size();
}
int g_ParserIntSequence(std::string str, std::vector<int>& vect)
{
    std::stringstream ss(str);
    int i;

    while (ss >> i)
    {
        vect.push_back(i);
        if (ss.peek() == ';')
            ss.ignore();
    }

    return vect.size();
}

// definitions of CCSVParser member functions
void CCSVParser::ConvertLineStringValueToIntegers()
{
    LineIntegerVector.clear();
    for (unsigned i = 0; i < LineFieldsValue.size(); ++i)
    {
        string si = LineFieldsValue[i];
        int value = atoi(si.c_str());

        if (value >= 1)
            LineIntegerVector.push_back(value);
    }
}

bool CCSVParser::OpenCSVFile(string fileName, bool b_required)
{
    mFileName = fileName;
    inFile.open(fileName.c_str());

    if (inFile.is_open())
    {
        if (IsFirstLineHeader)
        {
            string s;
            std::getline(inFile, s);
            vector<string> FieldNames = ParseLine(s);

            for (size_t i = 0;i < FieldNames.size();i++)
            {
                string tmp_str = FieldNames.at(i);
                size_t start = tmp_str.find_first_not_of(" ");

                string name;
                if (start == string::npos)
                {
                    name = "";
                }
                else
                {
                    name = tmp_str.substr(start);
                    //TRACE("%s,", name.c_str());
                }
                FieldsIndices[name] = (int)i;
            }
        }
        return true;
    }
    else
    {
        if (b_required)
        {
            dtalog.output() << "File " << fileName << " does not exist. Please check." << std::endl;
            //g_program_stop();
        }
        return false;
    }
}

bool CCSVParser::ReadRecord()
{
    LineFieldsValue.clear();

    if (inFile.is_open())
    {
        string s;
        std::getline(inFile, s);
        if (s.length() > 0)
        {
            LineFieldsValue = ParseLine(s);
            return true;
        }
        else
        {
            return false;
        }
    }
    else
    {
        return false;
    }
}

bool CCSVParser::ReadSectionHeader(string s)
{
    //skip // data
    Headers.clear();
    FieldsIndices.clear();

    if (s.length() == 0)
        return true;

    vector<string> FieldNames = ParseLine(s);

    for (size_t i = 0; i < FieldNames.size(); i++)
    {
        string tmp_str = FieldNames.at(i);
        size_t start = tmp_str.find_first_not_of(" ");

        string name;
        if (start == string::npos)
        {
            name = "";
        }
        else
        {
            name = tmp_str.substr(start);
        }
        Headers.push_back(name);
        FieldsIndices[name] = (int)i;
    }
    return true;
}

bool CCSVParser::ReadRecord_Section()
{
    LineFieldsValue.clear();

    if (inFile.is_open())
    {
        string s;
        std::getline(inFile, s);
        if (s.length() > 0)
        {
            if(s.find("[") != string::npos)  // synchro single csv file
            {
                LineFieldsValue = ParseLine(s);

                if (LineFieldsValue.size() >= 1)
                {
                    SectionName = LineFieldsValue[0];
                }

                //re-read section header
                ReadSectionHeader(s);
                std::getline(inFile, s);
            }
            LineFieldsValue = ParseLine(s);
            return true;
        }
        else
        {
            if (m_bLastSectionRead)  // reach the last section
                return false;
            else
            {
                if (inFile.eof())
                    return false;
                else
                    return true;
            }
        }
    }
    else
    {
        return false;
    }
}

vector<string> CCSVParser::ParseLine(string line)
{
    vector<string> SeperatedStrings;
    string subStr;

    if (line.length() == 0)
        return SeperatedStrings;

    std::istringstream ss(line);

    if (line.find_first_of('"') == string::npos)
    {
        while (std::getline(ss, subStr, Delimiter))
        {
            SeperatedStrings.push_back(subStr);
        }

        if (line.at(line.length() - 1) == ',')
        {
            SeperatedStrings.push_back("");
        }
    }
    else
    {
        while (line.length() > 0)
        {
            size_t n1 = line.find_first_of(',');
            size_t n2 = line.find_first_of('"');

            if (n1 == string::npos && n2 == string::npos) //last field without double quotes
            {
                subStr = line;
                SeperatedStrings.push_back(subStr);
                break;
            }

            if (n1 == string::npos && n2 != string::npos) //last field with double quotes
            {
                size_t n3 = line.find_first_of('"', n2 + 1); // second double quote

                //extract content from double quotes
                subStr = line.substr(n2 + 1, n3 - n2 - 1);
                SeperatedStrings.push_back(subStr);

                break;
            }

            if (n1 != string::npos && (n1 < n2 || n2 == string::npos))
            {
                subStr = line.substr(0, n1);
                SeperatedStrings.push_back(subStr);
                if (n1 < line.length() - 1)
                {
                    line = line.substr(n1 + 1);
                }
                else // comma is the last char in the line string, push an empty string to the back of vector
                {
                    SeperatedStrings.push_back("");
                    break;
                }
            }

            if (n1 != string::npos && n2 != string::npos && n2 < n1)
            {
                size_t n3 = line.find_first_of('"', n2 + 1); // second double quote
                subStr = line.substr(n2 + 1, n3 - n2 - 1);
                SeperatedStrings.push_back(subStr);
                size_t idx = line.find_first_of(',', n3 + 1);

                if (idx != string::npos)
                {
                    line = line.substr(idx + 1);
                }
                else
                {
                    break;
                }
            }
        }
    }
    return SeperatedStrings;
}

bool CCSVParser::GetValueByFieldName(string field_name, string& value, bool required_field)
{
    if (FieldsIndices.find(field_name) == FieldsIndices.end())
    {
        if (required_field)
        {
            dtalog.output() << "Field " << field_name << " in file " << mFileName << " does not exist. Please check the file." << std::endl;
            g_program_stop();
        }
        return false;
    }
    else
    {
        if (LineFieldsValue.size() == 0)
        {
            return false;
        }

        unsigned int index = FieldsIndices[field_name];
        if (index >= LineFieldsValue.size())
        {
            return false;
        }
        string str_value = LineFieldsValue[index];

        if (str_value.length() <= 0)
        {
            return false;
        }

        value = str_value;
        return true;
    }
}

bool g_read_a_line(FILE* f)
/* read a line from the current line from the file */
{

    char ch;

    while (1) {
        ch = getc(f);
        if (ch != 13 && ch != 10 && ch != EOF)
        {
            // do nothing
        }
        else { /* terminate if it's end of line or end of file */
            {
                // do nothing
            }
            if (ch == EOF)
                return false;

            return true;
        }
    }
}

double g_calculate_p2p_distance_in_meter_from_latitude_longitude(double longitud1, double latitud1, double longitud2, double latitud2)
{
    double PI = 3.14159265358979323846;
    double RADIO_TERRESTRE = 6372797.56085;
    double GRADOS_RADIANES = PI / 180;

    double haversine;
    double temp;
    double distancia_puntos;

    latitud1 = latitud1 * GRADOS_RADIANES;
    longitud1 = longitud1 * GRADOS_RADIANES;
    latitud2 = latitud2 * GRADOS_RADIANES;
    longitud2 = longitud2 * GRADOS_RADIANES;

    haversine = (pow(sin((1.0 / 2) * (latitud2 - latitud1)), 2)) + ((cos(latitud1)) * (cos(latitud2)) * (pow(sin((1.0 / 2) * (longitud2 - longitud1)), 2)));
    temp = 2 * asin(fmin(1.0, sqrt(haversine)));
    distancia_puntos = RADIO_TERRESTRE * temp;

    return distancia_puntos;
}




/* ***************************************************************************** */
/* Copyright:      Francois Panneton and Pierre L'Ecuyer, University of Montreal */
/*                 Makoto Matsumoto, Hiroshima University                        */
/* Notice:         This code can be used freely for personal, academic,          */
/*                 or non-commercial purposes. For commercial purposes,          */
/*                 please contact P. L'Ecuyer at: lecuyer@iro.UMontreal.ca       */
/* ***************************************************************************** */
#define W 32
#define R 16
#define P 0
#define M1 13
#define M2 9
#define M3 5

#define MAT0POS(t,v) (v^(v>>t))
#define MAT0NEG(t,v) (v^(v<<(-(t))))
#define MAT3NEG(t,v) (v<<(-(t)))
#define MAT4NEG(t,b,v) (v ^ ((v<<(-(t))) & b))

#define V0            STATE[state_i                   ]
#define VM1           STATE[(state_i+M1) & 0x0000000fU]
#define VM2           STATE[(state_i+M2) & 0x0000000fU]
#define VM3           STATE[(state_i+M3) & 0x0000000fU]
#define VRm1          STATE[(state_i+15) & 0x0000000fU]
#define VRm2          STATE[(state_i+14) & 0x0000000fU]
#define newV0         STATE[(state_i+15) & 0x0000000fU]
#define newV1         STATE[state_i                 ]
#define newVRm1       STATE[(state_i+14) & 0x0000000fU]

#define FACT 2.32830643653869628906e-10

static unsigned int state_i = 0;
static unsigned int STATE[R];
static unsigned int z0, z1, z2;
unsigned int g_RandomSeed = 100;
void InitWELLRNG512a(unsigned int* init) {
    int j;
    state_i = 0;
    for (j = 0; j < R; j++)
        STATE[j] = init[j];
}

double WELLRNG512a(void) {
    z0 = VRm1;
    z1 = MAT0NEG(-16, V0) ^ MAT0NEG(-15, VM1);
    z2 = MAT0POS(11, VM2);
    newV1 = z1 ^ z2;
    newV0 = MAT0NEG(-2, z0) ^ MAT0NEG(-18, z1) ^ MAT3NEG(-28, z2) ^ MAT4NEG(-5, 0xda442d24U, newV1);
    state_i = (state_i + 15) & 0x0000000fU;
    return ((double)STATE[state_i]) * FACT;
}

double g_get_random_ratio()
{
    //	g_RandomSeed = (g_LCG_a * g_RandomSeed + g_LCG_c) % g_LCG_M;  //m_RandomSeed is automatically updated.
    //	return float(g_RandomSeed)/g_LCG_M;

    return WELLRNG512a();
}

void generate_default_settings()
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
	myfile << "[assignment] ,,assignment_mode,column_generation_iterations,column_updating_iterations,odme_iterations,sensitivity_analysis_iterations,simulation_iterations,remarks" << '\n';
	myfile << ",,dta,5,0,0,0,1," << '\n';
	myfile << "[agent_type],agent_type_no,agent_type,name,display_code,vot,flow_type,pce,person_occupancy,desired_speed_ratio,headway,real_time_info," << '\n';
	myfile << ",1,auto,auto,auto,10,0,1,1,1,1,0," << '\n';
	myfile << ",2,walk,walk,walk,10,0,1,1,0.1,1,0," << '\n';
	myfile << ",3,bike,bike,bike,10,0,,1,0.2,1,0" << '\n';
	myfile << ",4,bus,bus,bus,10,0,1,10,0.75,5,0" << '\n';
	myfile << ",5,truck,truck,truck,10,0,1,1,1,3,0" << '\n';
	myfile << ",6,cav,auto,auto,10,0,1,1,1,1,0," << '\n';
	myfile << ",7,ev,ev,ev,10,0,1,1,1,1,0," << '\n';
	myfile << ",8,traveler,traveler,traveler,10,0,1,1,1,1,0," << '\n';
	myfile << "" << '\n';
	myfile << "[link_type],link_type,link_type_name,,agent_type_blocklist,type_code,traffic_flow_code,vdf_type," << '\n';
	myfile << ",1,motorway,,,f,0,qvdf," << '\n';
	myfile << ",2,trunk,,,a,0,qvdf," << '\n';
	myfile << ",3,primary,,,a,0,bpr," << '\n';
	myfile << ",4,residential,,,a,0,bpr," << '\n';
	myfile << ",5,secondary,,,a,0,bpr," << '\n';
	myfile << ",6,tertiary,,,a,0,bpr," << '\n';
	myfile << ",7,unclassified,,,a,0,bpr," << '\n';
	myfile << "[demand_period],demand_period_id,demand_period,,time_period" << '\n';
	myfile << ", 1,AM, ,0900_0930,, , , " << '\n';
	myfile << "[departure_time_profile], departure_time_profile_no, , ,time_period, , , ,T0000,T0005,T0010,T0015,T0020,T0025,T0030,T0035,T0040,T0045,T0050,T0055,T0100,T0105,T0110,T0115,T0120,T0125,T0130,T0135,T0140,T0145,T0150,T0155,T0200,T0205,T0210,T0215,T0220,T0224,T0230,T0235,T0239,T0245,T0250,T0255,T0300,T0305,T0310,T0315,T0320,T0325,T0330,T0335,T0340,T0345,T0350,T0355,T0400,T0404,T0410,T0415,T0420,T0425,T0430,T0435,T0440,T0445,T0450,T0455,T0500,T0505,T0510,T0515,T0520,T0525,T0530,T0535,T0540,T0545,T0550,T0555,T0600,T0605,T0610,T0615,T0620,T0625,T0630,T0635,T0640,T0645,T0650,T0655,T0700,T0705,T0710,T0715,T0720,T0725,T0730,T0735,T0740,T0745,T0750,T0755,T0800,T0805,T0810,T0815,T0820,T0825,T0830,T0835,T0840,T0845,T0850,T0855,T0900,T0905,T0910,T0915,T0920,T0925,T0930,T0935,T0940,T0945,T0950,T0955,T1000,T1005,T1010,T1015,T1020,T1025,T1030,T1035,T1040,T1045,T1050,T1055,T1100,T1105,T1110,T1115,T1120,T1125,T1130,T1135,T1140,T1145,T1150,T1155,T1200,T1205,T1210,T1215,T1220,T1225,T1230,T1235,T1240,T1245,T1250,T1255,T1300,T1305,T1310,T1315,T1320,T1325,T1330,T1335,T1340,T1345,T1350,T1355,T1400,T1405,T1410,T1415,T1420,T1425,T1430,T1435,T1440,T1445,T1450,T1455,T1500,T1505,T1510,T1515,T1520,T1525,T1530,T1535,T1540,T1545,T1550,T1555,T1600,T1605,T1610,T1615,T1620,T1625,T1630,T1635,T1640,T1645,T1650,T1655,T1700,T1705,T1710,T1715,T1720,T1725,T1730,T1735,T1740,T1745,T1750,T1755,T1800,T1805,T1810,T1815,T1820,T1825,T1830,T1835,T1840,T1845,T1850,T1855,T1900,T1905,T1910,T1915,T2020,T1925,T1930,T1935,T1940,T1945,T1950,T1955,T2000,T2005,T2010,T2015,T2020,T2025,T2030,T2035,T2040,T2045,T2050,T2055,T2100,T2105,T2110,T2115,T2120,T2125,T2130,T2135,T2140,T2145,T2150,T2155,T2200,T2205,T2210,T2215,T2220,T2225,T2230,T2235,T2240,T2245,T2250,T2255,T2300,T2305,T2310,T2315,T2320,T2325,T2330,T2335,T2340,T2345,T2350,T2355,T2400" << '\n';
	myfile << ", 1, , ,0900_0930, , , ,0.000571,0.000571,0.000571,0.000571,0.000571,0.000571,0.000506,0.000506,0.000506,0.000445,0.000445,0.000445,0.000391,0.000391,0.000391,0.000357,0.000357,0.000357,0.000328,0.000328,0.000328,0.000319,0.000319,0.000319,0.000302,0.000302,0.000302,0.000292,0.000292,0.000292,0.000296,0.000296,0.000296,0.00031,0.00031,0.00031,0.00031,0.00031,0.00031,0.000319,0.000319,0.000319,0.000383,0.000383,0.000383,0.000496,0.000496,0.000496,0.000568,0.000568,0.000568,0.000656,0.000656,0.000656,0.00095,0.00095,0.00095,0.001368,0.001368,0.001368,0.001587,0.001587,0.001587,0.00175,0.00175,0.00175,0.002288,0.002288,0.002288,0.002921,0.002921,0.002921,0.003242,0.003242,0.003242,0.003218,0.003218,0.003218,0.003803,0.003803,0.003803,0.004459,0.004459,0.004459,0.005002,0.005002,0.005002,0.005207,0.005207,0.005207,0.005677,0.005677,0.005677,0.005994,0.005994,0.005994,0.006018,0.006018,0.006018,0.005508,0.005508,0.005508,0.00529,0.00529,0.00529,0.005058,0.005058,0.005058,0.004833,0.004833,0.004833,0.004421,0.004421,0.004421,0.004327,0.004327,0.004327,0.004364,0.004364,0.004364,0.004343,0.004343,0.004343,0.004139,0.004139,0.004139,0.004201,0.004201,0.004201,0.004291,0.004291,0.004291,0.00435,0.00435,0.00435,0.004409,0.004409,0.004409,0.004566,0.004566,0.004566,0.004674,0.004674,0.004674,0.004761,0.004761,0.004761,0.004827,0.004827,0.004827,0.004882,0.004882,0.004882,0.0049,0.0049,0.0049,0.004887,0.004887,0.004887,0.004835,0.004835,0.004835,0.004899,0.004899,0.004899,0.005023,0.005023,0.005023,0.005065,0.005065,0.005065,0.005162,0.005162,0.005162,0.005436,0.005436,0.005436,0.005772,0.005772,0.005772,0.005907,0.005907,0.005907,0.005877,0.005877,0.005877,0.00605,0.00605,0.00605,0.006196,0.006196,0.006196,0.006248,0.006248,0.006248,0.006308,0.006308,0.006308,0.006404,0.006404,0.006404,0.006391,0.006391,0.006391,0.006401,0.006401,0.006401,0.006526,0.006526,0.006526,0.006574,0.006574,0.006574,0.006271,0.006271,0.006271,0.005937,0.005937,0.005937,0.005578,0.005578,0.005578,0.005293,0.005293,0.005293,0.004834,0.004834,0.004834,0.004387,0.004387,0.004387,0.00403,0.00403,0.00403,0.003748,0.003748,0.003748,0.003382,0.003382,0.003382,0.003121,0.003121,0.003121,0.002963,0.002963,0.002963,0.00289,0.00289,0.00289,0.002671,0.002671,0.002671,0.002468,0.002468,0.002468,0.002365,0.002365,0.002365,0.002249,0.002249,0.002249,0.002015,0.002015,0.002015,0.001784,0.001784,0.001784,0.00164,0.00164,0.00164,0.001474,0.001474,0.001474,0.001312,0.001312,0.001312,0.001132,0.001132,0.001132,0.001005,0.001005,0.001005,0.000889,0.000889,0.000889,0.000778,0.000778,0.000778,0.000676" << '\n';
	myfile << "[demand_file_list],file_sequence_no,file_name,,format_type,demand_period,agent_type,scale_factor,departure_time_profile_no," << '\n';
	myfile << ",1,input_matrix.csv,,matrix,AM,auto,0.2,1," << '\n';
	myfile << ",2,input_matrix.csv,,matrix,AM,walk,0.2,1," << '\n';
	myfile << ",3,input_matrix.csv,,matrix,AM,bike,0.2,1," << '\n';
	myfile << ",4,input_matrix.csv,,matrix,AM,bus,0.2,1," << '\n';
	myfile << ",5,input_matrix.csv,,matrix,AM,truck,0.2,1," << '\n';
	myfile << ",6,input_matrix.csv,,matrix,AM,cav,0.2,1," << '\n';
	myfile << ",7,input_matrix.csv,,matrix,AM,ev,0.2,1," << '\n';
	myfile << ",8,input_demand_column.csv,,column,AM,auto,0,1," << '\n';
	myfile << ",9,input_demand_path.csv,,path,AM,auto,0,1," << '\n';
	//myfile << ",8,input_activity_plan.csv,,activity_plan,AM,traveler,0.2,1," << '\n';
	//myfile << ",9,input_activity_plan.csv,,activity_plan,AM,ev,0.2,1," << '\n';

	myfile.close();
}

// perform network assignment using settings from settings.csv
void perform_network_assignment()
{
    // reset all the log files to defult 0: not output; if want to output these logs set to 1
	dtalog.output() << "DTALite Log" << std::fixed << std::setw(12) << '\n';
	dtalog.debug_level() = 0;
	dtalog.log_sig() = 0;
	dtalog.log_odme() = 0;
	dtalog.log_path() = 0;
	dtalog.log_dta() = 0;
	dtalog.log_ue() = 0;

	int assignment_mode = 1;
    int column_generation_iterations = 20;
	int column_updating_iterations = 40;
	int ODME_iterations = 20;
	int sensitivity_analysis_iterations = 0;
    int simulation_iterations = 0;
	int number_of_memory_blocks = 8;

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

				// these are the assignment modes
				// two usually methods are ue (user equilibrium) and dta (dynamic traffic assignment)
				// the main difference of these two methods are different output in link_performance.csv
				// for basic uses set assignment mode to 'ue'
				// for more detailed link performances (one minute) set 'dta'
				if (assignment_mode_str == "lue")
					assignment_mode = 0;
				else if (assignment_mode_str == "dta")
					assignment_mode = 1;
				else if (assignment_mode_str == "zone2access")  //access link
					assignment_mode = 21;
				else if (assignment_mode_str == "cbi")  //congestion bottleneck identification
					assignment_mode = 11;
				else if (assignment_mode_str == "cbsa")  //congestion bottleneck sensitivity
					assignment_mode = 12;
				else
				{
					dtalog.output() << "assignment_mode " << assignment_mode_str.c_str() << " in settings.csv is invalid.\n";
					g_program_stop();
				}

				// iteration number of reassignment
				parser_settings.GetValueByFieldName("column_updating_iterations", column_updating_iterations, true, true);

				if (assignment_mode == 1)
				{
					parser_settings.GetValueByFieldName("odme_iterations", ODME_iterations, true, false);
					parser_settings.GetValueByFieldName("sensitivity_analysis_iterations", sensitivity_analysis_iterations, true, false);
				}
				parser_settings.GetValueByFieldName("simulation_iterations", simulation_iterations, true, false);

				// the start interation of generating signals, if there is no signals set this number larger than the iteration number
				parser_settings.GetValueByFieldName("number_of_memory_blocks", number_of_memory_blocks, false, false);

				dtalog.output() << "number_of_memory_blocks = " << number_of_memory_blocks << " in settings.csv.\n";

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

    if (assignment_mode == 21)
    {
        dtalog.output() << "assignment mode zone2access is not implemented.\n";
        g_program_stop;
    }

    if (assignment_mode == 11)
    {
        perform_cbi();
        return;
    }

    if (assignment_mode == 12)
    {
        perform_cbsa();
        return;
    }

    network_assignment(assignment_mode, column_generation_iterations, column_updating_iterations, ODME_iterations, 
                       sensitivity_analysis_iterations, simulation_iterations, number_of_memory_blocks);
}