#include "bit_if.h"


unsigned int find_entry(vector<vector<int> > unique_COO, vector<int> idx){
    unsigned int nnz_TVM = unique_COO.size();
    int d = idx.size();
    int i_key = 0; 

    for (int i = 0; i < nnz_TVM; i++) {
        bool match = true;
        // Check if the first d-1 entries in the row match with the entries in v
        for (int j = 0; j < d; j++) {
            //cout<< unique_COO[i][j] <<" "<<idx[j]<<endl;
            if (unique_COO[i][j] != idx[j]) {
                match = false;
                continue;
            }
        }
        // If a match is found, return the index of the row
        if (match) {
            i_key = i; 
            break;
        }
    }
    return i_key; 
}

BitIf::BitIf(COO &coo){

    vector<vector<int>> coo_indices = coo.getIndices();
    vector<double> coo_val = coo.getVal();

    this->nnz = coo_val.size();
    this->dim = coo_indices[0].size();
    this->mode = coo.getMode(); 

    this->delta_d = new int*[dim]();
    this->val = new double[nnz]();
    this->bit_encoder.resize(dim*nnz);
    unsigned int *incr_counter = new unsigned int[dim]();
    this->num_increments = new unsigned int[dim]();

    for(int i = 0; i<dim; i++){
        unsigned int count = 1; 
        for(unsigned int row = 0; row <nnz-1; row++){
            int delta_j = coo_indices[row+1][i] - coo_indices[row][i];
            if (delta_j !=0)
                count +=1;      
        }
        this->num_increments[i] = count;   
    }

    this->val[0] = coo_val[0];
    for (int j = 0; j < dim; j++){
        this->delta_d[j] = new int[num_increments[j]]();
        this->delta_d[j][0] = coo_indices[0][j];
        incr_counter[j] = 1;
        bit_encoder[j] = 1;
    }

    for(int row = 1; row<nnz; row++){
        this->val[row] = coo_val[row];
        for(int j = 0; j< dim ; j++){
            int delta_j = coo_indices[row][j] - coo_indices[row-1][j];
            if (delta_j !=0){
                this->delta_d[j][incr_counter[j]] = delta_j;
                incr_counter[j]++;
                this->bit_encoder[dim*row+j] = 1;
            }
        }
    }
}

void BitIf::printIncrementArrays(){

    cout<<"Increment arrays:"<<"\n"; 
    for(int i=0;i<this->dim;i++){
        cout<<"["; 

        for(int j=0;j< this->num_increments[i] ;j++ ){
            cout<<delta_d[i][j]<<" ";}
            cout<<"]"<<"\n";
    }
    cout<<"\n";
}

void BitIf::printBitEncoder(){
    cout<<"Bit Encoding array:"<<"\n"<< "["; 
    for(int i=0;i<this->bit_encoder.size();i++)
        cout<<this->bit_encoder[i]<<" ";   
    cout<<"]"<<"\n\n";
}

void BitIf::computeTVM(vector<double>& vec){

    std::chrono::time_point<std::chrono::high_resolution_clock> TVMstart, TVMend, TVMend2;

    unsigned int i_k = this->delta_d[this->mode-1][0]; 
    double temp_val = vec[i_k]*this->val[0]; 
    int sum_bitset = 0; 
    unsigned int count_k = 1;

    int dim_TVM = this->dim - 1;
    this->TVMdelta_d = new int*[dim_TVM];

    vector<unsigned int> ind_count(dim); 
    fill(ind_count.begin(), ind_count.end(), 1);
    vector<int> idx(dim_TVM); 

    unordered_map<unsigned int, double> TVM_tensor_map;
    TVM_tensor_map.reserve(nnz);
    
    vector<vector<int> > TVM_COO; 

    unsigned int i_key = 0; 
    
    vector<unsigned int> loop_indices;
    for (unsigned int j = 0; j < this->dim; j++) {
        if (j != this->mode-1) {
            loop_indices.push_back(j);
        }
    }
    TVMstart = std::chrono::high_resolution_clock::now();
    int j_ = 0; 
    for(auto j : loop_indices){
        idx[j_] = delta_d[j][0];
        j_+=1;
    }

    unsigned int precomputed_hash = hash_idx(idx); 
    //vector<unsigned int> precomputed_hashes;
    //precomputed_hashes.reserve(nnz);
    //precomputed_hashes.push_back(hash_idx(idx));

    for (unsigned int i = 1; i < nnz; i++) {
        sum_bitset = 0;
        if (this->bit_encoder[i * dim + (mode - 1)]) {
            i_k += this->delta_d[this->mode - 1][count_k];
            count_k++;
        }

        for (auto j : loop_indices) {
            if (this->bit_encoder[i * dim + j]) {
                sum_bitset++;
            }
        }

        if (sum_bitset >= 1) {
            i_key = precomputed_hash; //precomputed_hashes.back();
            auto mapIterator = TVM_tensor_map.find(i_key);
            if (mapIterator != TVM_tensor_map.end()) {
                mapIterator->second += temp_val;
            } else {
                TVM_tensor_map.emplace(i_key, temp_val);
                TVM_COO.push_back(idx);
            }
            temp_val = vec[i_k] * this->val[i];
            j_ = 0;
            for (auto j : loop_indices) {
                if (this->bit_encoder[i * dim + j]) {
                    idx[j_] += this->delta_d[j][ind_count[j]];
                    ind_count[j] += 1;
                }
                j_ += 1;
            }
            //precomputed_hashes.push_back(hash_idx(idx));
            precomputed_hash = hash_idx(idx); 
        } else {
            temp_val += vec[i_k] * this->val[i];
        }

    }

    i_key = precomputed_hash; //precomputed_hashes.back();

    auto mapIterator = TVM_tensor_map.find(i_key);
    if (mapIterator != TVM_tensor_map.end()) {
        mapIterator->second += temp_val;
    } else {
        TVM_tensor_map.emplace(i_key, temp_val);
        TVM_COO.push_back(idx);
    }

    TVMend = std::chrono::high_resolution_clock::now();
    cout<< chrono::duration_cast<chrono::milliseconds>(TVMend - TVMstart).count() << " "; 
    
    if (TVM_COO.size() != TVM_tensor_map.size()) {
            throw std::runtime_error(" New COO and TVM map are not the same size");  
    }  
    for(auto indices: TVM_COO){
        unsigned int key = hash_idx(indices);
        this->TVMval.push_back(TVM_tensor_map[key]);
    }     

    unsigned int nnzTVM = TVMval.size();
    int dimTVM = TVM_COO[0].size();
    this->TVMdelta_d = new int*[dimTVM]();
    this->TVMbit_encoder.resize(dimTVM*nnzTVM);
    unsigned int *TVMincr_counter = new unsigned int[dimTVM]();
    unsigned int *num_incrementsTVM = new unsigned int[dimTVM]();

    for(int i = 0; i<dimTVM; i++){
    unsigned int count = 1; 
        for(unsigned int row = 0; row <nnzTVM-1; row++){
            int delta_j = TVM_COO[row+1][i] - TVM_COO[row][i];
            if (delta_j !=0)
                count +=1;      
        }
        num_incrementsTVM[i] = count;   
    }

    for (int j = 0; j < dimTVM; j++){
        this->TVMdelta_d[j] = new int[num_incrementsTVM[j]]();
        this->TVMdelta_d[j][0] = TVM_COO[0][j];
        TVMincr_counter[j] = 1;
        TVMbit_encoder[j] = 1;
    }

    for(int row = 1; row < nnzTVM; row++){
        for(int j = 0; j< dimTVM ; j++){
            int delta_j = TVM_COO[row][j] - TVM_COO[row-1][j];
            if (delta_j !=0){
                this->TVMdelta_d[j][TVMincr_counter[j]] = delta_j;
                TVMincr_counter[j]++;
                this->TVMbit_encoder[dimTVM*row+j] = 1;
            }
        }
    }

    TVMend2 = std::chrono::high_resolution_clock::now();
    cout<< chrono::duration_cast<chrono::milliseconds>(TVMend2 - TVMstart).count() <<endl; 
} 

int BitIf::getDim(){
    return this->dim; 
}

unsigned int* BitIf::getNumIncrements(){
    return this->num_increments;
}

void BitIf::validate(vector<double> val_test){
    unsigned int nnz = val_test.size(); 
    int ind_k = 0; 
    if (nnz != this->TVMval.size()) {
        throw std::runtime_error("Test failed: both vectors must have the same size");  
    }    
    double norm2 = 0.0; 
    double checksum = 0.0; 
    double inf_norm = 0.0; 
    double curr_val_test = 0.0;
    double curr_tvm_val = 0.0; 
    unsigned int incorrect_count = 0; 

    for (unsigned int row = 0; row < nnz; row++) {
        curr_val_test = val_test[row];
        curr_tvm_val = this->TVMval[row];
        if(abs(curr_tvm_val-curr_val_test)>0.1){
            incorrect_count +=1;
        }
            
        norm2 += pow((curr_tvm_val-curr_val_test),2);
        checksum += abs(curr_tvm_val-curr_val_test);
        inf_norm = max(abs(curr_tvm_val-curr_val_test), inf_norm);

    }
    norm2 = sqrt(norm2);
    cout << "2-norm: "<< norm2<< " "<< "checksum: "<< checksum<<" " << "#incorrect: "<< incorrect_count<<endl;
    cout << "2-norm: "<< norm2<< " "<< "checksum: "<< checksum<<" " << "inf norm: " << inf_norm << " "<<"#incorrect: "<< incorrect_count<<endl;
}

BitIf::~BitIf(){

    for (int i = 0; i < this->dim; i++){
        if(this->delta_d[i] != nullptr)
            delete [] this->delta_d[i];
    } 
    if(this->delta_d != nullptr)   
        delete [] this->delta_d;
    
    if (this->TVMdelta_d != nullptr){
        for (int i = 0; i < this->dim-1; i++) 
            delete [] this->TVMdelta_d[i];
        delete [] this->TVMdelta_d;
    }
    if (this->val != nullptr)
        delete [] this->val;
    if (this->num_increments != nullptr)
        delete [] this->num_increments;

    this->delta_d = nullptr;
    this->TVMdelta_d = nullptr;
    this->num_increments = nullptr;
    this->val = nullptr;
}

BlockedBitIf::BlockedBitIf(BlockedCOO &blocked_coo){

    vector<vector<vector<int> > > blocked_indices = blocked_coo.getBlockedIndices(); 
    vector<vector<double> > blocked_vals = blocked_coo.getBlockedVal();
    this->num_blocks = blocked_indices.size(); 
    this->dim = blocked_coo.getDim();
    this->mode = blocked_coo.getMode();
    this->n_d = blocked_coo.getTensorSizes();
    this->block_size = blocked_coo.getBlockSize();

    this->size_per_block = new unsigned int[num_blocks]();
    this->delta_d = new int**[num_blocks]();
    this->val = new double*[num_blocks]();
    this->num_increments = new unsigned int*[num_blocks]();
    this->bit_encoder.resize(num_blocks);
    unsigned int total_entries = 0; 
    unsigned int *total_incr = new unsigned int[dim];

    for(int i = 0; i<dim; ++i)
        total_incr[i]=0;

    //Allocate memory for blocked increment arrays and blocked tensor value arrays
    for (unsigned int b = 0; b < num_blocks; b++) {
        size_per_block[b] = blocked_indices[b].size();
        unsigned int size_b = size_per_block[b];
        total_entries+=size_b;
        if (size_b == 0)
            continue;       

        delta_d[b] = new int*[dim]();
        val[b] = new double[size_b]();
        unsigned int *incr_counter_b = new unsigned int[dim]();
        num_increments[b] = new unsigned int[dim]();
        bit_encoder[b].resize(dim*size_b);

        //compute the number of increments
        for(int i = 0; i<dim; i++){
            unsigned int count = 1; 
            for(unsigned int row = 0; row <size_b-1; row++){
                int delta_j = blocked_indices[b][row+1][i] - blocked_indices[b][row][i];
                if (delta_j !=0)
                    count +=1;      
                } 
            num_increments[b][i] = count;
            total_incr[i]+=count;
        }

        //allocate memory for the increments and set the inital values for delta_d
        this->val[b][0] = blocked_vals[b][0];
        unsigned int delta_d_size; 
        for (int j = 0; j < dim; j++){
            delta_d_size = num_increments[b][j];
            delta_d[b][j] = new int[delta_d_size]();
            delta_d[b][j][0] = blocked_indices[b][0][j];
            incr_counter_b[j] = 1;
            bit_encoder[b][j] = 1;
        }
        
        //add the increments and bit encodings in delta_d for each block
        for(int row = 1; row<size_b; row++){
            this->val[b][row] = blocked_vals[b][row];
            for(int j = 0; j< dim ; j++){
                int delta_j = blocked_indices[b][row][j] - blocked_indices[b][row-1][j];
                if (delta_j !=0){
                    delta_d[b][j][incr_counter_b[j]] = delta_j;
                    incr_counter_b[j]++;
                    bit_encoder[b][dim*row+j] = 1;
                }
            }
        }
        if(incr_counter_b !=nullptr)
            delete [] incr_counter_b;   
    }
}

int BlockedBitIf::getDim(){
    return this->dim;
}

void BlockedBitIf::computeTVMBlockedBytes(vector<double> &vec){
     
    int sum_bitset; 
    unsigned int i_key; 

    int TVMdim = dim-1; 
    vector<unsigned int> ind_count(dim); 
    vector<int> idx(TVMdim); 
    vector<unsigned int> loop_indices;
    for (unsigned int j = 0; j < this->dim; j++) {
        if (j != this->mode-1) {
            loop_indices.push_back(j);
        }
    }

    this->TVMn_d.resize(TVMdim);
    int i_= 0; 
    for (auto i : loop_indices){
        TVMn_d[i_] = n_d[i];
        i_+=1;
    }

    unordered_map<unsigned int, unordered_map<unsigned int, double> > blocked_TVM_tensor_map;
    vector<vector<vector<int> > > TVM_COO;
    unordered_map<unsigned int, unsigned int> block_idx_map;

    double temp_val = 0;
    unsigned int count_k = 1;
    unsigned int i_k = 0;
    unsigned int block_i = 0; 

    unsigned int block_sizeTVM = block_size/2; 

    TVM_total_bytes = 0.0;
    for (unsigned int b = 0; b < this->num_blocks; b++) {
        int j_ = 0; 
        TVM_total_bytes += sizeof(int);
        for(auto j : loop_indices){
            idx[j_] = delta_d[b][j][0];
            j_ += 1;
            TVM_total_bytes += 2*sizeof(int); 
        }
        
        i_k = this->delta_d[b][this->mode-1][0];
        TVM_total_bytes += 2*sizeof(int); 
        temp_val = vec[i_k]*this->val[b][0];
        TVM_total_bytes += 3*sizeof(double);
        unsigned int size_b = this->size_per_block[b];
        TVM_total_bytes += 2*sizeof(int);

        if(size_b == 1){
            block_i = getBlock(idx, block_sizeTVM, TVMn_d);
            TVM_total_bytes += 2*sizeof(int) + TVMdim*sizeof(int);
            i_key = hash_idx(idx);
            TVM_total_bytes += sizeof(int) + TVMdim*sizeof(int);
            if (blocked_TVM_tensor_map.count(block_i) == 0) {
                blocked_TVM_tensor_map[block_i] = unordered_map<unsigned int, double>();
                TVM_total_bytes += sizeof(unordered_map<unsigned int, unordered_map<unsigned int, double> >) + sizeof(unordered_map<unsigned int, double>);
                unsigned int block_idx = TVM_COO.size();
                TVM_total_bytes += sizeof(int);
                TVM_COO.push_back(vector<vector<int>>());
                TVM_total_bytes += sizeof(vector<vector<int>>);
                block_idx_map[block_i] = block_idx;
                TVM_total_bytes += 2*sizeof(int);
            }
            unordered_map<unsigned int, double>& TVM_tensor_map = blocked_TVM_tensor_map[block_i];
            TVM_total_bytes += sizeof(int)*TVM_tensor_map.size() + sizeof(double)*TVM_tensor_map.size();
            
            if (TVM_tensor_map.count(i_key) == 0) {
                TVM_tensor_map[i_key] = temp_val;
                TVM_total_bytes += sizeof(double) + sizeof(int);
                unsigned int block_idx = block_idx_map[block_i];
                TVM_total_bytes += 2*sizeof(int);
                TVM_COO[block_idx].push_back(idx);
                TVM_total_bytes += sizeof(int) + sizeof(int)*TVMdim;
            } else {
                TVM_tensor_map[i_key] += temp_val;
                TVM_total_bytes += sizeof(double) + sizeof(int);
            }
            continue;
        }

        else{
            count_k = 1;
            fill(ind_count.begin(), ind_count.end(), 1);
            TVM_total_bytes += sizeof(int)*TVMdim;

            for (unsigned int i = 1; i< size_b; i++){
                sum_bitset = 0; 
                TVM_total_bytes += sizeof(int);
                if (this->bit_encoder[b][i*dim+(this->mode-1)]){
                    i_k+= this->delta_d[b][this->mode-1][count_k];
                    count_k++;     
                    TVM_total_bytes += 3*sizeof(int);
                }
                j_ = 0; 
                for (auto j : loop_indices) {
                    TVM_total_bytes += 1/8; //for accessing a bit
                    if (this->bit_encoder[b][i*dim+j]){
                        sum_bitset++;
                        TVM_total_bytes += sizeof(int);
                    }
                    j_+=1;      
                }   
                if (sum_bitset > 0){
                    block_i = getBlock(idx, block_sizeTVM, TVMn_d);
                    TVM_total_bytes += 2*sizeof(int) + TVMdim*sizeof(int);
                    i_key = hash_idx(idx);
                    TVM_total_bytes += sizeof(int) + TVMdim*sizeof(int);
                    if (blocked_TVM_tensor_map.count(block_i) == 0) {
                        blocked_TVM_tensor_map[block_i] = unordered_map<unsigned int, double>();
                        TVM_total_bytes += sizeof(unordered_map<unsigned int, unordered_map<unsigned int, double> >) + sizeof(unordered_map<unsigned int, double>);
                        unsigned int block_idx = TVM_COO.size();
                        TVM_total_bytes += sizeof(int);
                        TVM_COO.push_back(vector<vector<int>>());
                        TVM_total_bytes += sizeof(vector<vector<int>>);
                        block_idx_map[block_i] = block_idx;
                        TVM_total_bytes += 2*sizeof(int);
                    }
                    unordered_map<unsigned int, double>& TVM_tensor_map = blocked_TVM_tensor_map[block_i];
                    TVM_total_bytes += sizeof(int)*TVM_tensor_map.size() + sizeof(double)*TVM_tensor_map.size();
                    
                    if (TVM_tensor_map.count(i_key) == 0) {
                        TVM_tensor_map[i_key] = temp_val;
                        TVM_total_bytes += sizeof(double) + sizeof(int);
                        unsigned int block_idx = block_idx_map[block_i];
                        TVM_total_bytes += 2*sizeof(int);
                        TVM_COO[block_idx].push_back(idx);
                        TVM_total_bytes += sizeof(int) + sizeof(int)*TVMdim;
                    } else {
                        TVM_tensor_map[i_key] += temp_val;
                        TVM_total_bytes += sizeof(double) + sizeof(int);
                    }

                    temp_val = vec[i_k]*this->val[b][i];
                    TVM_total_bytes += 3*sizeof(double);
                    j_ = 0;
                    for (auto j : loop_indices) {
                        if (this->bit_encoder[b][i*dim+j]){
                            idx[j_] += this->delta_d[b][j][ind_count[j]];
                            ind_count[j] += 1;
                            TVM_total_bytes += 3*sizeof(int);
                        }
                        j_+=1; 
                    }
                }
                else{
                    temp_val += vec[i_k]*this->val[b][i];  
                    TVM_total_bytes += 3*sizeof(double);
                }  
            }
            block_i = getBlock(idx, block_sizeTVM, TVMn_d);
            TVM_total_bytes += 2*sizeof(int) + TVMdim*sizeof(int);
            i_key = hash_idx(idx);
            TVM_total_bytes += sizeof(int) + TVMdim*sizeof(int);
            if (blocked_TVM_tensor_map.count(block_i) == 0) {
                blocked_TVM_tensor_map[block_i] = unordered_map<unsigned int, double>();
                TVM_total_bytes += sizeof(unordered_map<unsigned int, unordered_map<unsigned int, double> >) + sizeof(unordered_map<unsigned int, double>);
                unsigned int block_idx = TVM_COO.size();
                TVM_total_bytes += sizeof(int);
                TVM_COO.push_back(vector<vector<int>>());
                //TVM_total_bytes += sizeof(vector<vector<int>>);
                block_idx_map[block_i] = block_idx;
                TVM_total_bytes += 2*sizeof(int);
            }
            unordered_map<unsigned int, double>& TVM_tensor_map = blocked_TVM_tensor_map[block_i];
            TVM_total_bytes += sizeof(int)*TVM_tensor_map.size() + sizeof(double)*TVM_tensor_map.size();

            if (TVM_tensor_map.count(i_key) == 0) {
                TVM_tensor_map[i_key] = temp_val;
                TVM_total_bytes += sizeof(double) + sizeof(int);
                unsigned int block_idx = block_idx_map[block_i];
                TVM_total_bytes += 2*sizeof(int);
                TVM_COO[block_idx].push_back(idx);
                TVM_total_bytes += sizeof(int) + sizeof(int)*TVMdim;
            } else {
                TVM_tensor_map[i_key] += temp_val;
                TVM_total_bytes += sizeof(double) + sizeof(int);
            }
        } 
    }

    if (TVM_COO.size() != blocked_TVM_tensor_map.size()) {
        throw std::runtime_error(" New COO and TVM map are not the same size");  
    }  
    this->TVMnum_blocks = TVM_COO.size(); 
    TVMval_blocked.resize(TVMnum_blocks);

    int beta_max = TVM_COO[0].size(), beta_min = TVM_COO[0].size();

    for (const auto& entry : blocked_TVM_tensor_map) {
        unsigned int block_i = entry.first;
        TVM_total_bytes += 2*sizeof(int);
        const unordered_map<unsigned int, double>& TVM_tensor_map = entry.second;
        TVM_total_bytes += sizeof(int)*TVM_tensor_map.size() + sizeof(double)*TVM_tensor_map.size();
        unsigned int block_idx = block_idx_map[block_i]; 
        TVM_total_bytes += 2*sizeof(int);
        int tvm_coo_size = TVM_COO[block_idx].size();
        beta_max = std::max(beta_max, tvm_coo_size);
        beta_min = std::min(beta_min, tvm_coo_size);
        for (auto indices : TVM_COO[block_idx]){
            unsigned int key = hash_idx(indices);
            TVM_total_bytes += sizeof(int) + TVMdim*sizeof(int);
            TVMval_blocked[block_idx].push_back(TVM_tensor_map.at(key));
            TVM_total_bytes += 2*sizeof(double);
        }
    }
    //convert it back to Bit-If structure:
    this-> TVMdim = this->dim - 1; 
    TVM_total_bytes += 2*sizeof(int);
    unsigned int** TVMnum_incr = new unsigned int*[TVMnum_blocks]();
    TVM_total_bytes += TVMnum_blocks*sizeof(int);
    this->TVMdelta_d_blocked = new int**[TVMnum_blocks]();
    TVM_total_bytes += TVMnum_blocks*sizeof(int);
    this->TVMbit_encoder_blocked.resize(TVMnum_blocks);
    TVM_total_bytes += TVMnum_blocks/8;


    for (unsigned int b = 0; b<TVMnum_blocks; ++b){
        unsigned int *incr_counter = new unsigned int[TVMdim]();
        //TVM_total_bytes += TVMdim*sizeof(int);
        unsigned int TVMsize_b = TVM_COO[b].size();
        TVM_total_bytes += 2*sizeof(int);
        TVMnum_incr[b] = new unsigned int[TVMdim]();
        TVM_total_bytes += TVMdim*sizeof(int);
        TVMdelta_d_blocked[b] = new int*[TVMdim]();
        TVM_total_bytes += TVMdim*sizeof(int);
        TVMbit_encoder_blocked[b].resize(TVMdim*TVMsize_b);
        TVM_total_bytes += TVMdim*TVMsize_b/8;

        for(int i = 0; i<TVMdim; i++){
            unsigned int count = 1; 
            for(unsigned int row = 0; row <TVMsize_b-1; row++){
                int delta_j = TVM_COO[b][row+1][i] - TVM_COO[b][row][i];
                TVM_total_bytes += sizeof(int);
                if (delta_j !=0){
                    count +=1;   
                    TVM_total_bytes += sizeof(int);
                }   
            } 
            TVMnum_incr[b][i] = count;
            TVM_total_bytes += sizeof(int);
        }

        for (int j = 0; j < TVMdim; j++){
            this->TVMdelta_d_blocked[b][j] = new int[TVMnum_incr[b][j]]();
            TVM_total_bytes += TVMnum_incr[b][j]*sizeof(int);
            this->TVMdelta_d_blocked[b][j][0] = TVM_COO[b][0][j];
            incr_counter[j] = 1;
            TVM_total_bytes += 3*sizeof(int);
            this->TVMbit_encoder_blocked[b][j] = 1;
            TVM_total_bytes += 2/8;
        }

        for(int row = 1; row < TVMsize_b; row++){
            for(int j = 0; j< TVMdim ; j++){
                int delta_j = TVM_COO[b][row][j] - TVM_COO[b][row-1][j];
                TVM_total_bytes += 2*sizeof(int);
                if (delta_j !=0){
                    this->TVMdelta_d_blocked[b][j][incr_counter[j]] = delta_j;
                    incr_counter[j]++;
                    TVM_total_bytes += 3*sizeof(int);
                    this->TVMbit_encoder_blocked[b][TVMdim*row+j] = 1;
                    TVM_total_bytes += 2/8;
                }
            }
        }
        if(incr_counter !=nullptr)
            delete [] incr_counter;
        
    }
    double beta_max_bytes = beta_max*TVMdim*sizeof(int) + beta_max*sizeof(double);
    double beta_min_bytes = beta_min*TVMdim*sizeof(int) + beta_min*sizeof(double);
    cout<<"total bytes moved (GB/s):"<<TVM_total_bytes/pow(2, 30)<<"; "<<"nnz max b:"<<beta_max<<"; "<<"bytes max b:"<<beta_max_bytes<<"; "<<"nnz min b:"<<beta_min<<"; "<<"bytes min b:"<<beta_min_bytes<<endl; 

    if(TVMnum_incr !=nullptr){
        for (unsigned int b = 0; b < this->TVMnum_blocks; b++) {
            delete [] TVMnum_incr[b];
        }
        delete [] TVMnum_incr;
    }
}

void BlockedBitIf::computeTVMBlockedRuntime(vector<double> &vec){
    
    std::chrono::time_point<std::chrono::high_resolution_clock> TVMstart, TVMend, midstart, midend, TVMend2;

    int sum_bitset; 
    unsigned int i_key; 

    int TVMdim = dim-1; 
    vector<unsigned int> ind_count(dim); 
    vector<int> idx(TVMdim); 

    vector<unsigned int> loop_indices;
    for (unsigned int j = 0; j < this->dim; j++) {
        if (j != this->mode-1) {
            loop_indices.push_back(j);
        }
    }

    this->TVMn_d.resize(TVMdim);
    int i_= 0; 
    for (auto i : loop_indices){
        TVMn_d[i_] = n_d[i];
        i_+=1;
    }

    unordered_map<unsigned int, unordered_map<unsigned int, double> > blocked_TVM_tensor_map;
    vector<vector<vector<int> > > TVM_COO;
    unordered_map<unsigned int, unsigned int> block_idx_map;

    double temp_val = 0;
    unsigned int count_k = 1;
    unsigned int i_k = 0;
    unsigned int block_i = 0; 
    unsigned int precomputed_hash;
    
    TVMstart = std::chrono::high_resolution_clock::now();

    for (unsigned int b = 0; b < this->num_blocks; b++) {
        
        i_k = this->delta_d[b][this->mode-1][0];
        temp_val = vec[i_k]*this->val[b][0];
        unsigned int size_b = this->size_per_block[b];
        unsigned int block_sizeTVM = block_size/2; 

        int j_ = 0; 
        for(auto j : loop_indices){
            idx[j_] = delta_d[b][j][0];
            j_ += 1;
        }
        precomputed_hash = hash_idx(idx);
        
        if(size_b == 1){
            
            block_i = getBlock(idx, block_sizeTVM, TVMn_d);
            i_key = precomputed_hash;
            if (blocked_TVM_tensor_map.count(block_i) == 0) {
                blocked_TVM_tensor_map[block_i] = unordered_map<unsigned int, double>();
                unsigned int block_idx = TVM_COO.size();
                TVM_COO.push_back(vector<vector<int>>());
                block_idx_map[block_i] = block_idx;
            }
            unordered_map<unsigned int, double>& TVM_tensor_map = blocked_TVM_tensor_map[block_i];
            
            if (TVM_tensor_map.count(i_key) == 0) {
                TVM_tensor_map[i_key] = temp_val;
                unsigned int block_idx = block_idx_map[block_i];
                TVM_COO[block_idx].push_back(idx);
            } else {
                TVM_tensor_map[i_key] += temp_val;
            }
            continue;
        }

        else{
            count_k = 1;
            fill(ind_count.begin(), ind_count.end(), 1);

            for (unsigned int i = 1; i< size_b; i++){
                sum_bitset = 0; 
                if (this->bit_encoder[b][i*dim+(this->mode-1)]){
                    i_k+= this->delta_d[b][this->mode-1][count_k];
                    count_k++;     
                }
                for (auto j : loop_indices) {
                    if (this->bit_encoder[b][i*dim+j]){
                        sum_bitset++;
                    }
                }   
                if (sum_bitset > 0){
                    block_i = getBlock(idx, block_sizeTVM, TVMn_d);
                    i_key = precomputed_hash;
                    if (blocked_TVM_tensor_map.count(block_i) == 0) {
                        blocked_TVM_tensor_map[block_i] = unordered_map<unsigned int, double>();
                        unsigned int block_idx = TVM_COO.size();
                        TVM_COO.push_back(vector<vector<int>>());
                        block_idx_map[block_i] = block_idx;
                    }
                    unordered_map<unsigned int, double>& TVM_tensor_map = blocked_TVM_tensor_map[block_i];
                    
                    if (TVM_tensor_map.count(i_key) == 0) {
                        TVM_tensor_map[i_key] = temp_val;
                        unsigned int block_idx = block_idx_map[block_i];
                        TVM_COO[block_idx].push_back(idx);
                    } else {
                        TVM_tensor_map[i_key] += temp_val;
                    }
                    temp_val = vec[i_k]*this->val[b][i];

                    j_ = 0; 
                    for (auto j : loop_indices) {
                        if (this->bit_encoder[b][i*dim+j]){
                            idx[j_] += this->delta_d[b][j][ind_count[j]];
                            ind_count[j] += 1;
                        }
                        j_+=1;      
                    }
                    precomputed_hash = hash_idx(idx);
                }
                else{
                    temp_val += vec[i_k]*this->val[b][i];  
                }  
            }
            block_i = getBlock(idx, block_sizeTVM, TVMn_d);
            i_key = precomputed_hash;
            if (blocked_TVM_tensor_map.count(block_i) == 0) {
                blocked_TVM_tensor_map[block_i] = unordered_map<unsigned int, double>();
                unsigned int block_idx = TVM_COO.size();
                TVM_COO.push_back(vector<vector<int>>());
                block_idx_map[block_i] = block_idx;
            }
            unordered_map<unsigned int, double>& TVM_tensor_map = blocked_TVM_tensor_map[block_i];

            if (TVM_tensor_map.count(i_key) == 0) {
                TVM_tensor_map[i_key] = temp_val;
                unsigned int block_idx = block_idx_map[block_i];
                TVM_COO[block_idx].push_back(idx);
            } else {
                TVM_tensor_map[i_key] += temp_val;
            }
        }
    }

    TVMend = std::chrono::high_resolution_clock::now();
    this->TVMruntime = chrono::duration_cast<chrono::milliseconds>(TVMend - TVMstart).count();
    cout<<this->TVMruntime <<" "; 

    if (TVM_COO.size() != blocked_TVM_tensor_map.size()) {
        throw std::runtime_error(" New COO and TVM map are not the same size");  
    }  
    this->TVMnum_blocks = TVM_COO.size(); 
    TVMval_blocked.resize(TVMnum_blocks);

    for (const auto& entry : blocked_TVM_tensor_map) {
        const unsigned int& block_i = entry.first;
        const unordered_map<unsigned int, double>& TVM_tensor_map = entry.second;
        const unsigned int& block_idx = block_idx_map[block_i]; 
        for (const auto& indices : TVM_COO[block_idx]){
            unsigned int key = hash_idx(indices);
            const auto &it = TVM_tensor_map.find(key);
            assert( it != TVM_tensor_map.end() );
            TVMval_blocked[block_idx].push_back(it->second);
        }
    }
    //convert it back to Bit-If structure:
    this-> TVMdim = this->dim - 1; 
    unsigned int** TVMnum_incr = new unsigned int*[TVMnum_blocks]();
    this->TVMdelta_d_blocked = new int**[TVMnum_blocks]();
    this->TVMbit_encoder_blocked.resize(TVMnum_blocks);

    for (unsigned int b = 0; b<TVMnum_blocks; ++b){
        unsigned int *incr_counter = new unsigned int[TVMdim]();
        unsigned int TVMsize_b = TVM_COO[b].size();
        TVMnum_incr[b] = new unsigned int[TVMdim]();
        TVMdelta_d_blocked[b] = new int*[TVMdim]();
        TVMbit_encoder_blocked[b].resize(TVMdim*TVMsize_b);

        for(int i = 0; i<TVMdim; i++){
            unsigned int count = 1; 
            for(unsigned int row = 0; row <TVMsize_b-1; row++){
                int delta_j = TVM_COO[b][row+1][i] - TVM_COO[b][row][i];
                if (delta_j !=0){
                    count +=1;   
                }   
            } 
            TVMnum_incr[b][i] = count;
        }

        for (int j = 0; j < TVMdim; j++){
            this->TVMdelta_d_blocked[b][j] = new int[TVMnum_incr[b][j]]();
            this->TVMdelta_d_blocked[b][j][0] = TVM_COO[b][0][j];
            incr_counter[j] = 1;
            this->TVMbit_encoder_blocked[b][j] = 1;
        }

        for(int row = 1; row < TVMsize_b; row++){
            for(int j = 0; j< TVMdim ; j++){
                int delta_j = TVM_COO[b][row][j] - TVM_COO[b][row-1][j];
                if (delta_j !=0){
                    this->TVMdelta_d_blocked[b][j][incr_counter[j]] = delta_j;
                    incr_counter[j]++;
                    this->TVMbit_encoder_blocked[b][TVMdim*row+j] = 1;
                }
            }
        }
        if(incr_counter !=nullptr)
            delete [] incr_counter;
    }
    if(TVMnum_incr !=nullptr){
        for (unsigned int b = 0; b < this->TVMnum_blocks; b++) {
            delete [] TVMnum_incr[b];
        }
        delete [] TVMnum_incr;
    }
    TVMend2 = std::chrono::high_resolution_clock::now();
    cout<<chrono::duration_cast<chrono::milliseconds>(TVMend2 - TVMstart).count() <<endl; 
    
}

void BlockedBitIf::computeTVMBytes(vector<double> &vec){
    TVM_total_bytes = 0.0; 

    int sum_bitset; 
    unsigned int size_in_b;
    unsigned int i_key; 

    int TVMdim = dim-1; 
    vector<unsigned int> ind_count(dim); 
    vector<int> idx(TVMdim); 
    vector<int> prev_idx(TVMdim);

    unordered_map<unsigned int, double> TVM_tensor_map;
    vector<vector<int> > TVM_COO; 

    bool new_entry = false; 
    double temp_val = 0;
    unsigned int count_k = 1;
    unsigned int i_k = 0;

    for (unsigned int b = 0; b < this->num_blocks; b++) {
        int j_ = 0; 
        TVM_total_bytes += sizeof(int);
        for(int j = 0; j<dim; ++j){
            if(j == this->mode-1)
                continue; 
            idx[j_] = delta_d[b][j][0];
            prev_idx[j_] = idx[j_];
            j_ += 1;
            TVM_total_bytes += 3*sizeof(int); 
        }
        
        i_k = this->delta_d[b][this->mode-1][0];
        TVM_total_bytes += 2*sizeof(int); 
        temp_val = vec[i_k]*this->val[b][0];
        TVM_total_bytes += 3*sizeof(double);
        size_in_b = this->size_per_block[b];
        TVM_total_bytes += 2*sizeof(int);

        if(size_in_b == 1){
            i_key = hash_idx(idx);
            TVM_total_bytes += sizeof(int) + TVMdim*sizeof(int);
            auto mapIterator = TVM_tensor_map.find(i_key);
            TVM_total_bytes += sizeof(int) + 2*sizeof(double);
            if (mapIterator != TVM_tensor_map.end()) {
                mapIterator->second += temp_val;
                TVM_total_bytes += sizeof(double) + sizeof(int);
            } else {
                TVM_tensor_map.insert({i_key, temp_val});
                TVM_total_bytes += sizeof(int) + 2*sizeof(double);//overhead is 32 
                TVM_COO.push_back(idx);
                TVM_total_bytes += sizeof(int)*TVMdim;
            } 
            continue;
        }

        else{
            count_k = 1;
            fill(ind_count.begin(), ind_count.end(), 1);
            TVM_total_bytes += sizeof(int)*(TVMdim+1);
            for (unsigned int i = 1; i< size_in_b; i++){
                sum_bitset = 0; 
                TVM_total_bytes += sizeof(int);
                if (this->bit_encoder[b][i*dim+(this->mode-1)]){
                    i_k+= this->delta_d[b][this->mode-1][count_k];
                    count_k++;     
                    TVM_total_bytes += 3*sizeof(int);
                }
                j_ = 0; 
                for (unsigned int j = 0; j < dim; j++) {
                    if(j == this->mode-1)
                        continue; 
                    TVM_total_bytes += 1/8;
                    if (this->bit_encoder[b][i*dim+j]){
                        sum_bitset++;
                        idx[j_] += this->delta_d[b][j][ind_count[j]];
                        ind_count[j] += 1;
                        TVM_total_bytes += 3*sizeof(int);
                    }
                    j_+=1;      
                }   
                if (sum_bitset > 0){
                    i_key = hash_idx(prev_idx);
                    TVM_total_bytes += sizeof(int) + TVMdim*sizeof(int);
                    auto mapIterator = TVM_tensor_map.find(i_key);
                    TVM_total_bytes += sizeof(int) + sizeof(double);
                    if (mapIterator != TVM_tensor_map.end()) {
                        mapIterator->second += temp_val;
                        TVM_total_bytes += sizeof(double) + sizeof(int);
                    } else {
                        TVM_tensor_map.insert({i_key, temp_val});
                        TVM_total_bytes += sizeof(int) + 2*sizeof(double);//overhead is 32 
                        TVM_COO.push_back(prev_idx);
                        TVM_total_bytes += sizeof(int)*TVMdim;
                    }     
                    temp_val = vec[i_k]*this->val[b][i];
                    TVM_total_bytes += 3*sizeof(double);
                    for (unsigned int j = 0; j < TVMdim; j++) {
                       prev_idx[j] = idx[j];
                       TVM_total_bytes += 2*sizeof(int);
                    }
                }
                else{
                    temp_val += vec[i_k]*this->val[b][i];  
                    TVM_total_bytes += 3*sizeof(double);
                }  
            }
            i_key = hash_idx(idx);
            TVM_total_bytes += sizeof(int) + TVMdim*sizeof(int);
            auto mapIterator = TVM_tensor_map.find(i_key);
            TVM_total_bytes += sizeof(int) + sizeof(double); //iterator has approx 8 bytes
            if (mapIterator != TVM_tensor_map.end()) {
                mapIterator->second += temp_val;
                TVM_total_bytes += sizeof(double) + sizeof(int);
            } else {
                TVM_tensor_map.insert({i_key, temp_val});
                TVM_total_bytes += sizeof(int) + 2*sizeof(double);//overhead is 32 
                TVM_COO.push_back(idx);
                TVM_total_bytes += sizeof(int)*TVMdim;
            } 
        }     
    }

    if (TVM_COO.size() != TVM_tensor_map.size()) {
        throw std::runtime_error(" New COO and TVM map are not the same size");  
    }  

    for(auto indices: TVM_COO){
        unsigned int key = hash_idx(indices);
        TVM_total_bytes += sizeof(int) + TVMdim*sizeof(int);
        TVMval.push_back(TVM_tensor_map[key]);
        TVM_total_bytes += 2*sizeof(double);
    }

    //convert it back to Bit-If structure:
    unsigned int nnz_TVM = TVMval.size();
    this-> TVMdim = this->dim - 1; 
    TVM_total_bytes += 4*sizeof(int);
    unsigned int* TVMnum_incr = new unsigned int[TVMdim](); 
    TVM_total_bytes += TVMdim*sizeof(int);
    this->TVMdelta_d = new int*[TVMdim]();
    TVM_total_bytes += TVMdim*sizeof(int);
    this->TVMbit_encoder.resize(TVMdim*nnz_TVM);
    TVM_total_bytes += TVMdim*nnz_TVM/8;
    unsigned int *incr_counter = new unsigned int[TVMdim]();
    TVM_total_bytes += TVMdim*sizeof(int);

    for(int i = 0; i<TVMdim; i++){
        unsigned int count = 1; 
        TVM_total_bytes += sizeof(int);
        for(unsigned int row = 0; row <nnz_TVM-1; row++){
            int delta_j = TVM_COO[row+1][i] - TVM_COO[row][i];
            TVM_total_bytes += 3*sizeof(int);
            if (delta_j !=0){
               count +=1; 
               TVM_total_bytes += sizeof(int);   
            }  
        } 
        TVMnum_incr[i] = count;   
        TVM_total_bytes += sizeof(int);
    }
   
    for (int j = 0; j < TVMdim; j++){
        this->TVMdelta_d[j] = new int[TVMnum_incr[j]]();
        TVM_total_bytes += TVMnum_incr[j]*sizeof(int);
        this->TVMdelta_d[j][0] = TVM_COO[0][j];
        incr_counter[j] = 1;
        TVM_total_bytes += 3*sizeof(int);
        this->TVMbit_encoder[j] = 1;
        TVM_total_bytes += 2/8;
    }

    for(int row = 1; row < nnz_TVM; row++){
        for(int j = 0; j< TVMdim ; j++){
            int delta_j = TVM_COO[row][j] - TVM_COO[row-1][j];
            TVM_total_bytes += 3*sizeof(int);
            if (delta_j !=0){
                this->TVMdelta_d[j][incr_counter[j]] = delta_j;
                incr_counter[j]++;
                TVM_total_bytes += 3*sizeof(int);
                this->TVMbit_encoder[TVMdim*row+j] = 1;
                TVM_total_bytes += 2/8;
            }
        }
    }
    if(incr_counter !=nullptr)
            delete [] incr_counter;
    
    cout<<TVM_total_bytes/pow(2, 30)<<endl; 
}

void BlockedBitIf::computeTVMRuntime(vector<double> &vec) {

    std::chrono::time_point<std::chrono::high_resolution_clock> TVMstart, TVMend;

    int sum_bitset; 
    unsigned int size_in_b;
    unsigned int i_key; 

    int TVMdim = dim-1; 
    vector<unsigned int> ind_count(dim); 
    vector<int> idx(TVMdim); 

    unordered_map<unsigned int, double> TVM_tensor_map;
    vector<vector<int> > TVM_COO; 

    bool new_entry = false; 
    double temp_val = 0;
    unsigned int count_k = 1;
    unsigned int i_k = 0;
    unsigned int precomputed_hash; 

    vector<unsigned int> loop_indices;
    for (unsigned int j = 0; j < this->dim; j++) {
        if (j != this->mode-1) {
            loop_indices.push_back(j);
        }
    }

    TVMstart = std::chrono::high_resolution_clock::now();

    for (unsigned int b = 0; b < this->num_blocks; b++) {
        int j_ = 0; 
        for(auto j : loop_indices){
            idx[j_] = delta_d[b][j][0];
            j_ += 1;
        }
        
        i_k = this->delta_d[b][this->mode-1][0];
        temp_val = vec[i_k]*this->val[b][0];
        size_in_b = this->size_per_block[b];
        precomputed_hash = hash_idx(idx);

        if(size_in_b == 1){
            i_key = precomputed_hash;
            auto mapIterator = TVM_tensor_map.find(i_key);
            if (mapIterator != TVM_tensor_map.end()) {
                mapIterator->second += temp_val;
            } else {
                TVM_tensor_map.emplace(i_key, temp_val);
                TVM_COO.push_back(idx);
            } 
            continue;
        }

        else{
            count_k = 1;
            fill(ind_count.begin(), ind_count.end(), 1);
            for (unsigned int i = 1; i< size_in_b; i++){
                sum_bitset = 0; 

                if (this->bit_encoder[b][i*dim+(this->mode-1)]){
                    i_k+= this->delta_d[b][this->mode-1][count_k];
                    count_k++;     
                } 
                for (auto j : loop_indices) {
                    if (this->bit_encoder[b][i*dim+j]){
                        sum_bitset++;
                    }    
                }   
                if (sum_bitset > 0){
                    i_key = precomputed_hash;
                    auto mapIterator = TVM_tensor_map.find(i_key);
                    if (mapIterator != TVM_tensor_map.end()) {
                        mapIterator->second += temp_val;
                    } else {
                        TVM_tensor_map.emplace(i_key, temp_val);
                        TVM_COO.push_back(idx);
                    }     
                    temp_val = vec[i_k]*this->val[b][i];

                    j_ = 0; 
                    for (auto j : loop_indices) {
                        if (this->bit_encoder[b][i*dim+j]){
                            idx[j_] += this->delta_d[b][j][ind_count[j]];
                            ind_count[j] += 1;
                        }
                        j_+=1;      
                    }
                    precomputed_hash = hash_idx(idx);
                }
                else{
                    temp_val += vec[i_k]*this->val[b][i];  
                }  
            }
            i_key = precomputed_hash;
            auto mapIterator = TVM_tensor_map.find(i_key);
            if (mapIterator != TVM_tensor_map.end()) {
                mapIterator->second += temp_val;
            } else {
                TVM_tensor_map.insert({i_key, temp_val});
                TVM_COO.push_back(idx);
            } 
        }
    }

    if (TVM_COO.size() != TVM_tensor_map.size()) {
        throw std::runtime_error(" New COO and TVM map are not the same size");  
    }  

    for(auto indices: TVM_COO){
        unsigned int key = hash_idx(indices);
        TVMval.push_back(TVM_tensor_map[key]);
    }

    //convert it back to Bit-If structure:
    unsigned int nnz_TVM = TVMval.size();
    this-> TVMdim = this->dim - 1; 
    unsigned int* TVMnum_incr = new unsigned int[TVMdim](); 
    this->TVMdelta_d = new int*[TVMdim]();
    this->TVMbit_encoder.resize(TVMdim*nnz_TVM);
    unsigned int *incr_counter = new unsigned int[TVMdim]();
    for(int i = 0; i<TVMdim; i++){
        unsigned int count = 1; 
        for(unsigned int row = 0; row <nnz_TVM-1; row++){
            int delta_j = TVM_COO[row+1][i] - TVM_COO[row][i];
            if (delta_j !=0)
               count +=1;      
        } 
        TVMnum_incr[i] = count;   
    }
   
    for (int j = 0; j < TVMdim; j++){
        this->TVMdelta_d[j] = new int[TVMnum_incr[j]]();
        this->TVMdelta_d[j][0] = TVM_COO[0][j];
        incr_counter[j] = 1;
        this->TVMbit_encoder[j] = 1;
    }

    for(int row = 1; row < nnz_TVM; row++){
        for(int j = 0; j< TVMdim ; j++){
            int delta_j = TVM_COO[row][j] - TVM_COO[row-1][j];
            if (delta_j !=0){
                this->TVMdelta_d[j][incr_counter[j]] = delta_j;
                incr_counter[j]++;
                this->TVMbit_encoder[TVMdim*row+j] = 1;
            }
        }
    }

    TVMend = std::chrono::high_resolution_clock::now();
    this->TVMruntime = chrono::duration_cast<chrono::milliseconds>(TVMend - TVMstart).count();
    cout<<this->TVMruntime<<endl;

}

void BlockedBitIf::validate(vector<double> val_test){

    unsigned int nnz = val_test.size(); 
    int ind_k = 0; 
    if (nnz != this->TVMval.size()) {
        throw std::runtime_error("Test failed: TVMval has different size than val_test");  
    }    
    double norm2 = 0.0; 
    double checksum = 0.0; 
    double inf_norm = 0.0; 
    double curr_val_test = 0.0;
    double curr_tvm_val = 0.0; 
    unsigned int incorrect_count = 0; 

    for (unsigned int row = 0; row < nnz; row++) {
        curr_val_test = val_test[row];
        curr_tvm_val = this->TVMval[row];
        if(abs(curr_tvm_val-curr_val_test)>0.1){
            incorrect_count+=1; 
        }
        norm2 += pow((curr_tvm_val-curr_val_test),2);
        checksum += abs(curr_tvm_val-curr_val_test);
        inf_norm = max(abs(curr_tvm_val-curr_val_test), inf_norm);
    }
    norm2 = sqrt(norm2);
    cout << "2-norm: "<< norm2<< " "<< "checksum: "<< checksum<<" " << "inf norm: " << inf_norm << " "<<"#incorrect: "<< incorrect_count<<endl;
}

void BlockedBitIf::validateBlocked(vector<vector<double> >val_test){

    unsigned int nnz = val_test.size(); 
    int ind_k = 0; 

    if (nnz != this->TVMval_blocked.size()) {
        throw std::runtime_error("Test failed: TVMval has different block size than val_test");  
    }    
    double norm2 = 0.0; 
    double checksum = 0.0; 
    double inf_norm = 0.0; 
    double curr_val_test = 0.0;
    double curr_tvm_val = 0.0; 
    unsigned int incorrect_count = 0; 

    for (unsigned int b = 0; b < TVMnum_blocks; b++){
        if (val_test[b].size() != this->TVMval_blocked[b].size()) {
            throw std::runtime_error("Test failed: TVMval has different size than val_test within block");  
        }  
        unsigned int TVMsize_b = TVMval_blocked[b].size();
        for (unsigned int row = 0; row < TVMsize_b; row++) {
            curr_val_test = val_test[b][row];
            curr_tvm_val = this->TVMval_blocked[b][row];
            if(abs(curr_tvm_val-curr_val_test)>0.1){
                incorrect_count+=1; 
            }
            norm2 += pow((curr_tvm_val-curr_val_test),2);
            checksum += abs(curr_tvm_val-curr_val_test);
            inf_norm = max(abs(curr_tvm_val-curr_val_test), inf_norm);
        }
    
    }
    norm2 = sqrt(norm2);
    cout << "2-norm: "<< norm2<< " "<< "checksum: "<< checksum<<" " << "inf norm: " << inf_norm << " "<<"#incorrect: "<< incorrect_count<<endl;
}

vector<vector<double>> BlockedBitIf::getValTVMBlocked()
{
    return this->getValTVMBlocked();
}

BlockedBitIf::~BlockedBitIf(){
    if(this->delta_d !=nullptr) {
        for (unsigned int b = 0; b < this->num_blocks; b++) {
            for(int j = 0; j<this->dim; j++){
                    delete [] this->delta_d[b][j];
            }
                delete [] this->delta_d[b];
        }
        delete [] this->delta_d;
    } 
    
    if(this->val !=nullptr){
        for (unsigned int b = 0; b < this->num_blocks; b++) {
            delete [] this->val[b];
        }
        delete [] this->val;
    } 

    //data after TVM is computed
    if(this->TVMdelta_d_blocked !=nullptr) {
        for (unsigned int b = 0; b < this->TVMnum_blocks; b++) {
            for(int j = 0; j<this->TVMdim; j++){
                    delete [] this->TVMdelta_d_blocked[b][j];
            }
                delete [] this->TVMdelta_d_blocked[b];
        }
        delete [] this->TVMdelta_d_blocked;
    } 

    if (this->TVMdelta_d != nullptr){ 
        for(int j = 0; j< this->TVMdim; j++){
            if(this->TVMdelta_d[j] != nullptr)
                delete [] this->TVMdelta_d[j];
        }
        delete [] this->TVMdelta_d; 
    }

    if(this->num_increments !=nullptr){
        for (unsigned int b = 0; b < this->num_blocks; b++) {
            delete [] this->num_increments[b];
        }
        delete [] this->num_increments;
    } 

    this->num_increments = nullptr;
    this->delta_d = nullptr; 
    this->val = nullptr;
    this->TVMdelta_d = nullptr; 
    this->TVMdelta_d_blocked = nullptr; 
}

vector<double> BlockedBitIf::getValTVM(){
    return this->TVMval;
}
