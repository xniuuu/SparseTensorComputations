#ifndef BITIF_H
#define BITIF_H

#include <vector>
#include <cstdint>
#include <cmath>
#include <algorithm>
#include <boost/dynamic_bitset.hpp>
#include "coo.h"

using namespace std;

class BitIf {
  public:
    BitIf(COO &coo);

    void printIncrementArrays();
    void printBitEncoder();
    void computeTVM(vector<double>& vec);

    int getDim(); 
    unsigned int* getNumIncrements();
    void validate(vector<double> val_test);

    //void removeBitIf();
    ~BitIf();
    
  private:
    unsigned int nnz;  
    int dim; 
    int mode; 
    int **delta_d = nullptr; //maybe replace int with int_fast32_t?
    double *val = nullptr;
    boost::dynamic_bitset<> bit_encoder;
    
    unsigned int *num_increments = nullptr;
    //vector<vector<int> > unique_COO; 
    
    //members to comute TVM:
    int** TVMdelta_d = nullptr;
    vector<double>  TVMval;
    boost::dynamic_bitset<> TVMbit_encoder;
};

class BlockedBitIf {

  public:
    BlockedBitIf(BlockedCOO &blocked_coo);

    //void computeTVM(double** vec);
    void computeTVMBlockedBytes(vector<double>& vec);
    void computeTVMBlockedRuntime(vector<double>& vec);

    void computeTVMBytes(vector<double>& vec);
    void computeTVMRuntime(vector<double>& vec);

    int getDim(); 
    void validate(vector<double>  val_test);
    void validateBlocked(vector<vector<double> > val_test);
    vector<vector<double> > getValTVMBlocked();
    vector<double> getValTVM();
    //void removeBitIf();
    ~BlockedBitIf();

  private:
    unsigned int block_size; 
    unsigned int TVMblock_size; 
    //mfor measuring bandwidth
    double TVM_total_bytes = 0.0; 
    double TVMruntime = 0.0; 
    
    int mode; 
    int ***delta_d = nullptr;
    double **val = nullptr;
    vector<boost::dynamic_bitset<>> bit_encoder;
    unsigned int num_blocks;
    unsigned int *size_per_block = nullptr; //stores the number of elements in each block length of array is num_blocks
    int dim;
    unsigned int **num_increments = nullptr;
    //vector<vector<int> > unique_COO; 
    vector<unsigned int> n_d; 

    //members to comute TVM:
    int **TVMdelta_d = nullptr;
    vector<double> TVMval;
    boost::dynamic_bitset<> TVMbit_encoder;

    vector<unsigned int> TVMn_d;
    //vector<vector<unsigned int> > TVMnum_incr;
    int ***TVMdelta_d_blocked = nullptr;
    vector<vector<double> > TVMval_blocked;
    vector<boost::dynamic_bitset<> > TVMbit_encoder_blocked;
    int TVMdim; 
    unsigned int TVMnum_blocks;
    
    
}; 

#endif
