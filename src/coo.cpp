#include "coo.h"

using namespace std; 

unsigned int hash_idx(const vector<int>& idx) {
  unsigned int seed = idx.size();
  for(auto x : idx) {
    x = ((x >> 16) ^ x) * 0x45d9f3b;
    x = ((x >> 16) ^ x) * 0x45d9f3b;
    x = (x >> 16) ^ x;
    seed ^= x + 0x9e3779b9 + (seed << 6) + (seed >> 2);
  }
  return seed;
}

unsigned int getBlock(const vector<int> &inds, unsigned int block_size,  vector<unsigned int> n_d){
    int num_dims = inds.size();
    std::vector<int> num_blocks_per_dim(num_dims);
    for (int i = 0; i < num_dims; ++i) {
        num_blocks_per_dim[i] = (n_d[i]%block_size == 0)? n_d[i]/block_size: n_d[i]/block_size + 1;
    }
    unsigned int block_number = 0;
    for (int i = 0; i < num_dims; ++i) {
        int block_index = inds[i] / block_size;
        int dim_multiplier = 1;
        for (int j = i + 1; j < num_dims; ++j) {
            dim_multiplier *= num_blocks_per_dim[j];
        }
        block_number += block_index * dim_multiplier;
    }
    return block_number;
}

unsigned int find_entry_coo(vector<vector<int> > unique_COO, vector<int> idx){
    unsigned int nnz_TVM = unique_COO.size();
    int d = idx.size();
    int i_key = 0; 

    for (int i = 0; i < nnz_TVM; i++) {
        bool match = true;
        for (int j = 0; j < d; j++) {
            if (unique_COO[i][j] != idx[j]) {
                match = false;
                continue;
            }
        }
        if (match) {
            i_key = i; 
            break;
        }
    }
    return i_key; 
}

unsigned int grayCode(unsigned int n) {
    return n ^ (n >> 1);
}

unsigned int inverseGrayCode(unsigned int n) {
    unsigned int i = 1;
    while (i < sizeof(unsigned int) * 8) {
        n ^= n >> i;
        i <<= 1;
    }
    return n;
}

COO::COO(const string file_name, bool offset, const char sep, vector<unsigned int> n_d, int mode){

    this->tensor_size = n_d;  
    this-> mode = mode;    
    this->dim = n_d.size();

    if (this->mode <= 0 || this->mode > dim) {
        throw std::runtime_error(" Invalid mode.");  
    }                     
    string line, elem;
    this->nnz = 0; 
    fstream file_tns (file_name, ios::in);
    if(file_tns.is_open()){        
        while(getline(file_tns, line)){
        this->nnz+=1;

        stringstream entry(line);
        vector<int> coo;
        vector<double> value;
            while(getline(entry, elem, sep)){
                coo.push_back(stoi(elem)-offset);
                value.push_back(stod(elem));
            }
        val.push_back(value.back());
        coo.pop_back(); // Remove the last element (column)
        indices.push_back(coo);
        assert(coo.begin() != coo.end());
        }
    }
    else
        throw std::runtime_error("Could not open the file\n");


}

unsigned int COO::returnZOrder(vector<int> idx){
    unsigned int z = 0;
    int max_n = 0;
    int d = idx.size();

    // Find the maximum index value across all dimensions
    for (int i = 0; i < d; i++) {
        max_n = max(max_n, idx[i]);
    }

    // Calculate the number of bits required to represent the maximum index value
    unsigned int num_bits = floor(log2(max_n)) + 1;

    // Iterate through each dimension and bit position
    // Interleave the bits of the indices to compute the z-order curve number
    for (int j = 0; j < num_bits; j++) {
        for (int i = d - 1; i >= 0; i--) {
            unsigned int bit = (idx[i] >> j) & 1;
            z |= bit << ((d-1-i) + d * j);
        }
    }

    return z;
}

unsigned int COO::returnHilbertNumber(vector<int> idx, int hilbert_curve_order){
    
    int d = idx.size();

    int max_n = 0;
    for (int i = 0; i < d; i++) {
        max_n = max(max_n, idx[i]);
    }

    unsigned int n = floor(log2(max_n));

       for (int i = 0; i < n; i++) {
        max_n = max(max_n, idx[i]);
    }

    unsigned int num_bits = floor(log2(max_n))+1;
    unsigned int hilbert_order = 0;

    unsigned int coord_bits = 0;

    for (int bit = num_bits - 1; bit >= 0; --bit) {
        unsigned int mask = 1 << bit;
        bool invert = false;
        unsigned int direction = 0;

        for (unsigned int coord : idx) {
            if (coord & mask) {
                direction ^= 1 << coord_bits;
                invert = !invert;
            }
        }

        if (!invert) {
            direction = ~direction;
        }

        unsigned int gray = grayCode(direction) >> (sizeof(unsigned int) * 8 - 1 - coord_bits);
        hilbert_order <<= coord_bits;
        hilbert_order |= gray;

        coord_bits += idx.size();
        while (coord_bits >= num_bits) {
            coord_bits -= num_bits;
            unsigned int tmp = hilbert_order;
            hilbert_order >>= coord_bits;
            hilbert_order ^= tmp;
            hilbert_order &= (1 << num_bits) - 1;
        }
    }

    return inverseGrayCode(hilbert_order);
}

vector<vector<int>> COO::getIndices()
{
    return indices;
}

void COO::printValues(){
    cout<<"Tensor values:"<<"\n"<< "["; 
    for(int i=0;i<val.size();i++)
        cout<<val[i]<<" ";
    std::cout<<"]"<<"\n";

}

vector<double> COO::getVal()
{
    return val;
}

vector<double> COO::getValTVM()
{
    return valTVM;
}

int COO::getMode()
{
    return this->mode;
}

void COO::computeTVM(vector<double> &vec){

    std::chrono::time_point<std::chrono::high_resolution_clock> TVMstart, TVMend, TVMend2;

    unordered_map<unsigned int, double> TVM_tensor_map;
    vector<vector<int> > TVM_COO; 

    int dim = this->indices[0].size();
    int TVMdim = dim-1;
    unsigned int i_k = this->indices[0][this->mode-1];
    double temp_val = this->val[0]*vec[i_k];
    int sum_incr = 0; 
    bool new_entry = false;
    unsigned int i_key = 0; 
    vector<int> idx(TVMdim);
    
    vector<unsigned int> loop_indices;
    for (unsigned int j = 0; j < this->dim; j++) {
        if (j != this->mode-1) {
            loop_indices.push_back(j);
        }
    }

    TVMstart = std::chrono::high_resolution_clock::now();

    for(unsigned int row = 1; row < this->nnz; row++){
        int i_ = 0; 
        for (int i = 0; i < dim; i++){  
            if(i == this->mode-1)
                continue;
            idx[i_] = indices[row-1][i];
            i_+=1; 
        }
        sum_incr = 0; 
        i_k = this->indices[row][this->mode-1]; 
        for (auto i : loop_indices){   
            if(this->indices[row][i] != this->indices[row-1][i]){
                sum_incr += 1; 
            }               
        }
        if (sum_incr > 0){
            i_key = hash_idx(idx);
            auto mapIterator = TVM_tensor_map.find(i_key);
            if (mapIterator != TVM_tensor_map.end()) {
                mapIterator->second += temp_val;
            } else {
                TVM_tensor_map.insert({i_key, temp_val});
                TVM_COO.push_back(idx);
            }
            temp_val = vec[i_k]*this->val[row];
        }else{
            temp_val += this->val[row]*vec[i_k];
        }             
    }
    int i_ = 0;
    for (auto i : loop_indices){  
        idx[i_] = indices[nnz-1][i];
        i_+=1; 
    }
    i_key = hash_idx(idx);
    auto mapIterator = TVM_tensor_map.find(i_key);
    if (mapIterator != TVM_tensor_map.end()) {
        mapIterator->second += temp_val;
    } else {
        TVM_tensor_map.insert({i_key, temp_val});
        TVM_COO.push_back(idx);
    }
    TVMend = std::chrono::high_resolution_clock::now();
    cout<<chrono::duration_cast<chrono::milliseconds>(TVMend - TVMstart).count()<<" ";
    //cout<<TVM_COO.size()<<" "<< TVM_tensor_map.size()<<endl;
    if (TVM_COO.size() != TVM_tensor_map.size()) {
        throw std::runtime_error(" New COO and TVM map are not the same size");  
    }  
    
    for(auto indices: TVM_COO){
        unsigned int key = hash_idx(indices);
        valTVM.push_back(TVM_tensor_map[key]);
    }
    TVMend2 = std::chrono::high_resolution_clock::now();
    cout<<chrono::duration_cast<chrono::milliseconds>(TVMend2 - TVMstart).count()<<endl;
}

/*-------------------------------------BLOCKED COO-----------------------------------------------*/

BlockedCOO::BlockedCOO(const string file_name_tns, bool offset, const char sep, vector<unsigned int> n_d, 
                        int b, int mode, const string traversal_curve){ 

    this->dim = n_d.size();
    this->tensor_size = n_d;     
    this->mode = mode;    
    if (this->mode <= 0 || this->mode > dim) {
        throw std::runtime_error(" Invalid mode. ");  
    }                
    
    this->block_size = b; 
    this->nnz = 0; 

    unsigned int num_blocks = 1;  
    unsigned int block_i = 0; 

    map<unsigned int, vector<vector<int> > > block_map_COO;
    map<unsigned int, vector<double > > block_map_val;

    vector<unsigned int> n_block(this->dim, 0);

    for (int i = 0; i < this->dim; i++) {
        unsigned int n_blocks_i = (this->tensor_size[i] + block_size - 1) / block_size;
        //cout<<this->tensor_size[i]<< " "<<  block_size<<endl; 
        num_blocks *= n_blocks_i;
        n_block[i] = n_blocks_i;
        //cout << n_blocks_i<<endl;
    }

    if (num_blocks > std::numeric_limits<int>::max()) {
        throw std::runtime_error(" The number of total blocks exceeds 2^{31}-1. Please choose a larger block size.");  
    }

    string line_tns, elem_tns;
    fstream file_tns (file_name_tns, ios::in);
    if(file_tns.is_open()){
        while(getline(file_tns, line_tns)){
            nnz +=1; 

            stringstream entry(line_tns);

            vector<int> coo;
            vector<double> value;
                while(getline(entry, elem_tns, sep)){
                    coo.push_back(stoi(elem_tns)-offset);
                    value.push_back(stod(elem_tns));
                }
            coo.pop_back(); // Remove the last element (column)
            block_i = getBlock(coo, this->block_size, this->tensor_size);
            block_map_COO[block_i].push_back(coo);
            block_map_val[block_i].push_back(value.back());
            
        }
    }
    else{
        throw std::runtime_error("Could not open the file\n");
    }

    if(traversal_curve == "hierarchical"){
        //this takes too lonk dang it, try find a better solution
        for (auto it = block_map_COO.begin(); it != block_map_COO.end(); ++it) {
            if (!it->second.empty()) {
                blocked_indices.push_back(it->second);
            }
        }
        for (auto it = block_map_val.begin(); it != block_map_val.end(); ++it) {
            if (!it->second.empty()) {
                blocked_val.push_back(it->second);
            }
        }

    } else if(traversal_curve == "zcurve"){
        sortBlockedCOOZCurve(n_block, block_map_COO, block_map_val);
    } else if(traversal_curve == "hilbert"){
        sortBlockedCOOHilbert(n_block, block_map_COO, block_map_val);
    } else{
        throw std::runtime_error("Invalid traversing curve");   
    }
}

unsigned int BlockedCOO::returnZOrder(vector<unsigned int> n_block, unsigned int block_number){
    
    unsigned int z = 0; 
    unsigned int max_n = 0;
    unsigned max_index = 0; 

    vector<unsigned int> block_indices(this->dim, 0);

    for (int i = 0; i < this->dim; i++){
        block_indices[i] = block_number%n_block[i];
        block_number/=n_block[i];
    }

    for (int i = 0; i < this->dim; i++){
        max_n = max(max_n, n_block[i]);
    }

    unsigned int  numBits = ceil(log2(max_n));
    for (int i = this->dim - 1; i >= 0; i--) {
        for (int j = numBits - 1; j >= 0; j--) {
            if ((block_indices[i] >> j) & 1) {
                z |= 1u << ((this->dim - 1 - i) * numBits + j);
            }
        }
    }
    return z; 
}

void BlockedCOO::sortBlockedCOOZCurve(vector<unsigned int> n_block, map<unsigned int, vector<vector<int>>> &block_map_COO, map<unsigned int, vector<double>> &block_map_val){

    bool are_keys_equal = true;
    for (const auto& kv : block_map_COO) {
        if (block_map_val.find(kv.first) == block_map_val.end()) {
            are_keys_equal = false;
            break;
        }
    }
    if (are_keys_equal) {
        for (const auto& kv : block_map_val) {
            if (block_map_COO.find(kv.first) == block_map_COO.end()) {
                are_keys_equal = false;
                break;
            }
        }
    }

    if (!are_keys_equal) {
        throw std::runtime_error("Blocked COO storage and associated blocked value storage does not match.");  
    }

    unsigned int nneb = block_map_COO.size(); //nneb = number of non-empy blocks
    blocked_indices.resize(nneb);
    blocked_val.resize(nneb);

    unsigned int z_number = 0; 
    map<unsigned int, vector<unsigned int> > z_order;
    
    for (const auto& block : block_map_val){
        unsigned int block_i = block.first;
        z_number = returnZOrder(n_block, block_i);
        if (z_order.find(z_number) == z_order.end()) {
            z_order[z_number] = vector<unsigned int>{block_i};
        } else {
            z_order[z_number].push_back(block_i);
        }
    }

    int total_z_entries = 0;
    for (const auto& entry : z_order) {
        total_z_entries += entry.second.size();
    }
    if (!(nneb == total_z_entries)) {
        throw std::runtime_error("Ordering map of z curve does not match coo block map");  
    }

    unsigned int count = 0; 
    for (const auto& entry : z_order){
        const vector<unsigned int>& block_indices = entry.second;
        for (const auto& block_i : block_indices) {
            blocked_indices_2d.insert(blocked_indices_2d.end(), block_map_COO[block_i].begin(), block_map_COO[block_i].end());
            blocked_val_2d.insert(blocked_val_2d.end(), block_map_val[block_i].begin(), block_map_val[block_i].end());

            blocked_indices[count] = block_map_COO[block_i];
            blocked_val[count] = block_map_val[block_i];
            count += 1;
        }
    }
}

unsigned int BlockedCOO::returnHilbertNumber(vector<unsigned int> n_block, unsigned int block_number, unsigned int hilbert_curve_order){

    unsigned int max_n = 0;
    vector<unsigned int> block_indices(this->dim, 0);

    for (int i = 0; i < this->dim; i++){
        block_indices[i] = block_number%n_block[i];
        block_number/=n_block[i];
    }

    int n = block_indices.size();
    
    for (int i = 0; i < n; i++) {
        max_n = max(max_n, block_indices[i]);
    }

    unsigned int num_bits = floor(log2(max_n))+1;
    unsigned int hilbert_order = 0;

    unsigned int coord_bits = 0;

    for (int bit = num_bits - 1; bit >= 0; --bit) {
        unsigned int mask = 1 << bit;
        bool invert = false;
        unsigned int direction = 0;

        for (unsigned int coord : block_indices) {
            if (coord & mask) {
                direction ^= 1 << coord_bits;
                invert = !invert;
            }
        }

        if (!invert) {
            direction = ~direction;
        }

        unsigned int gray = grayCode(direction) >> (sizeof(unsigned int) * 8 - 1 - coord_bits);
        hilbert_order <<= coord_bits;
        hilbert_order |= gray;

        coord_bits += block_indices.size();
        while (coord_bits >= num_bits) {
            coord_bits -= num_bits;
            unsigned int tmp = hilbert_order;
            hilbert_order >>= coord_bits;
            hilbert_order ^= tmp;
            hilbert_order &= (1 << num_bits) - 1;
        }
    }

    return inverseGrayCode(hilbert_order);
}

void BlockedCOO::sortBlockedCOOHilbert(vector<unsigned int> n_block, map<unsigned int, vector<vector<int>>> &block_map_COO, map<unsigned int, vector<double>> &block_map_val){

    bool are_keys_equal = true;
    for (const auto& kv : block_map_COO) {
        if (block_map_val.find(kv.first) == block_map_val.end()) {
            are_keys_equal = false;
            break;
        }
    }
    if (are_keys_equal) {
        for (const auto& kv : block_map_val) {
            if (block_map_COO.find(kv.first) == block_map_COO.end()) {
                are_keys_equal = false;
                break;
            }
        }
    }

    if (!are_keys_equal) {
        throw std::runtime_error("Blocked COO storage and associated blocked value storage does not match.");  
    }

    unsigned int nneb = block_map_COO.size(); //nneb = number of non-empy blocks
    
    blocked_indices.resize(nneb);
    blocked_val.resize(nneb);

    unsigned int hilbert_number = 0; 
    map<unsigned int, vector<unsigned int> > hilbert_order;

    unsigned int hilbert_curve_order = 1; 

    for (const auto& block : block_map_val){
        unsigned int block_i = block.first;
        hilbert_number = returnHilbertNumber(n_block, block_i, hilbert_curve_order);
    
        if (hilbert_order.find(hilbert_number) == hilbert_order.end()) {
            hilbert_order[hilbert_number] = vector<unsigned int>{block_i};
        } else {
            hilbert_order[hilbert_number].push_back(block_i);
        }
    }


    int total_hilbert_entries = 0;
    for (const auto& entry : hilbert_order) {
        total_hilbert_entries += entry.second.size();
    }
    if (!(nneb == total_hilbert_entries)) {
        throw std::runtime_error("Ordering map of z curve does not match coo block map");  
    }

    unsigned int count = 0; 
    for (const auto& num : hilbert_order){
        const vector<unsigned int>& block_indices = num.second;
        for (const auto& block_i : block_indices) {
            blocked_indices_2d.insert(blocked_indices_2d.end(), block_map_COO[block_i].begin(), block_map_COO[block_i].end());
            blocked_val_2d.insert(blocked_val_2d.end(), block_map_val[block_i].begin(), block_map_val[block_i].end());

            blocked_indices[count] = block_map_COO[block_i];
            blocked_val[count] = block_map_val[block_i];
            count += 1;
        }
    }
}

vector<vector<vector<int>>> BlockedCOO::getBlockedIndices(){
    return this->blocked_indices;
}

vector<unsigned int> BlockedCOO::getTensorSizes(){
    return this->tensor_size;
}

vector<vector<double>> BlockedCOO::getBlockedVal(){
    return this->blocked_val;
}

vector< vector<double > > BlockedCOO::getValTVM(){
    return valTVM;
}

unsigned int BlockedCOO::getBlockSize(){
    return this->block_size;
}

int BlockedCOO::getMode(){
    return this->mode;
}

int BlockedCOO::getDim(){
    return this->dim;
}

void BlockedCOO::computeTVM(vector<double> &vec){

    std::chrono::time_point<std::chrono::high_resolution_clock> TVMstart, TVMend, TVMend2;

    unsigned int num_blocks = blocked_val.size();
    int dim = this->blocked_indices[0][0].size();
    int TVMdim = dim-1; 

    unsigned int i_k = this->blocked_indices[0][0][this->mode-1];
    double temp_val = this->blocked_val[0][0]*vec[i_k];
    int sum_incr = 0; 
    unsigned int i_key = 0; 
    unsigned int block_i = 0; 
    unsigned int block_sizeTVM = block_size/2; 
    vector<int> idx(TVMdim);

    this->TVMtensor_size.resize(TVMdim);
    int i_= 0; 
    for (int i = 0; i < dim; i++){
        if(i == this->mode-1)
            continue;
        TVMtensor_size[i_] = tensor_size[i];
        i_+=1;
    }

    unordered_map<unsigned int, unordered_map<unsigned int, double> > blocked_TVM_tensor_map;
    vector<vector<vector<int> > > TVM_COO;
    unordered_map<unsigned int, unsigned int> block_idx_map;

    TVMstart = std::chrono::high_resolution_clock::now();

    for (unsigned int b = 0; b < num_blocks; b++) {
        i_k = this->blocked_indices[b][0][this->mode-1];
        temp_val = this->blocked_val[b][0]*vec[i_k];
        unsigned int size_b = this->blocked_indices[b].size();
        
        for(unsigned int row = 1; row < size_b; row++){
            int i_ = 0; 
            for (int i = 0; i < dim; i++){  
                if(i == this->mode-1)
                    continue;
                idx[i_] = blocked_indices[b][row-1][i];
                i_+=1; 
            }
            sum_incr = 0; 
            i_k = this->blocked_indices[b][row][this->mode-1];
            for (int i = 0; i < dim; i++){   
                if (i==this->mode-1)
                    continue;
                if(this->blocked_indices[b][row][i] != this->blocked_indices[b][row-1][i]){
                    sum_incr += 1; 
                }               
            }
            if (sum_incr > 0){
                i_key = hash_idx(idx);
                block_i = getBlock(idx, block_sizeTVM, TVMtensor_size);
                if (blocked_TVM_tensor_map.count(block_i) == 0) {
                    blocked_TVM_tensor_map[block_i] = unordered_map<unsigned int, double>();
                    unsigned int block_idx = TVM_COO.size();
                    TVM_COO.push_back(vector<vector<int> >());
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
                temp_val = vec[i_k]*this->blocked_val[b][row];
            }else{
                temp_val += this->blocked_val[b][row]*vec[i_k];
            }             
        }
        int i_ = 0; 
        for (int i = 0; i < dim; i++){  
            if(i == this->mode-1)
                continue;
            idx[i_] = blocked_indices[b][size_b-1][i];
            i_+=1; 
        }
        block_i = getBlock(idx, block_sizeTVM, TVMtensor_size);
        i_key = hash_idx(idx);
        if (blocked_TVM_tensor_map.count(block_i) == 0) {
            blocked_TVM_tensor_map[block_i] = unordered_map<unsigned int, double>();
            unsigned int block_idx = TVM_COO.size();
            TVM_COO.push_back(vector<vector<int> >());
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
    TVMend = std::chrono::high_resolution_clock::now();
    cout<<chrono::duration_cast<chrono::milliseconds>(TVMend - TVMstart).count()<<" ";
    if (TVM_COO.size() != blocked_TVM_tensor_map.size()) {
        throw std::runtime_error(" New COO and TVM map are not the same size");  
    }  
    unsigned int TVMnum_blocks = TVM_COO.size();
    valTVM.resize(TVMnum_blocks);

    for (const auto& entry : blocked_TVM_tensor_map) {
        unsigned int block_i = entry.first;
        const unordered_map<unsigned int, double>& TVM_tensor_map = entry.second;
        unsigned int block_idx = block_idx_map[block_i]; 
        for (auto indices : TVM_COO[block_idx]){
            unsigned int key = hash_idx(indices);
            valTVM[block_idx].push_back(TVM_tensor_map.at(key));
        }
    }
    
    TVMend2= std::chrono::high_resolution_clock::now();
    cout<<chrono::duration_cast<chrono::milliseconds>(TVMend2 - TVMstart).count()<<endl;

}
