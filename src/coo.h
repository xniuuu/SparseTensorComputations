#ifndef COO_H
#define COO_H

#include <cmath>
#include <vector>
#include <stdexcept>
#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include <sstream>
#include <algorithm>
#include <limits>
#include <cassert>
#include <unordered_map>
#include <chrono>

using namespace std;

//Hash function from: https://stackoverflow.com/questions/20511347/a-good-hash-function-for-a-vector/72073933#72073933
unsigned int hash_idx(const vector<int>& idx);
unsigned int getBlock(const vector<int> &inds, unsigned int block_size,  vector<unsigned int> n_d);

class COO {
  public:
    COO(const string file_name, bool offset, const char sep, vector<unsigned int> n_d, int mode);

    unsigned int returnZOrder(vector<int> idx);
    unsigned int returnHilbertNumber(vector<int> idx, int hilbert_curve_order);
    void computeTVM(vector<double>& vec);
    vector<vector<int> > getIndices(); 

    vector<double>  getVal();
    vector<vector<int> > getUniqueCOO(); 
    vector<double>  getValTVM();
    int getMode(); 
    void printValues();

  private:
    vector<vector<int>> indices;
    vector<double> val;
    vector<unsigned int> tensor_size; 
    int dim;
    unsigned int nnz; 
    int mode; 

    vector<double> valTVM;
};

class BlockedCOO {
  public:

    BlockedCOO(const string file_name_tns, bool offset, const char sep, vector<unsigned int> n_d, int b, int mode, const string traversal_curve);

    void sortBlockedCOOZCurve(vector<unsigned int> n_block, map<unsigned int, vector<vector<int> > >& block_map_COO, 
                            map<unsigned int, vector<double> >& block_map_val);
    unsigned int returnZOrder(vector<unsigned int> n_block, unsigned int block_number);

    void sortBlockedCOOHilbert(vector<unsigned int> n_block, map<unsigned int, vector<vector<int> > >& block_map_COO, 
                            map<unsigned int, vector<double> >& block_map_val);
    unsigned int returnHilbertNumber(vector<unsigned int> n_block, unsigned int block_number, unsigned int hilbert_curve_order);

    void computeTVM(vector<double> &vec);

    vector<vector<vector<int> > > getBlockedIndices(); 
    vector<unsigned int > getTensorSizes(); 
    vector<vector<double> > getBlockedVal();
    vector< vector<double > > getValTVM();
    unsigned int getBlockSize(); 
    int getMode(); 
    int getDim(); 

  private:
    vector<vector<vector<int> > > blocked_indices;
    vector<vector<double> > blocked_val;

    vector<vector<double> > valTVM;

    double** blocked_vec_ptr;
    vector<vector<int>> blocked_indices_2d;
    vector<double> blocked_val_2d;
    vector<unsigned int> tensor_size; 
    vector<unsigned int> TVMtensor_size; 

    unsigned int block_size; 
    int dim;
    int mode; 
    unsigned int nnz; 
};


#endif
