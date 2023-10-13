#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include <sstream>
#include "bit_if.h"
#include "coo.h"
//#include "block_map.h"

using namespace std;

//Using std::vector
vector<vector<double>> csv2COO(const string file_name, char sep);
//Using Armadillo
void print_arr(vector<vector<double>> content);

vector<vector<int> > csvToCOOIndices(const string file_name, char sep);
vector<double> csvToCOOValues(const string file_name, char sep);
vector<double> getVecFromFile(const string file_name);

unsigned int* determineIncrements(vector<vector<int>> COO); 
